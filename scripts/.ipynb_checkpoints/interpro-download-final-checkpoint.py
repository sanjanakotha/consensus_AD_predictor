#!/usr/bin/env python3

import sys, json, ssl, os, argparse
from urllib import request
from urllib.error import HTTPError
from time import sleep
from concurrent.futures import ThreadPoolExecutor, as_completed

HEADER_SEPARATOR = "|"
LINE_LENGTH = 80
MAX_WORKERS = 32  # Adjust based on system/network capacity

def fetch_sequences(ipr_id, reviewed=True):
    print(f"Starting {ipr_id}")
    status = "reviewed" if reviewed else "unreviewed"
    filename = f"../output/interpro/{status}/{ipr_id}.fasta"
    if os.path.exists(filename):
        print(f"[⏩] Skipped {ipr_id} (already exists at {filename})")
        return

    BASE_URL = f"https://www.ebi.ac.uk:443/interpro/api/protein/{status}/entry/InterPro/{ipr_id}/?page_size=200&extra_fields=sequence"
    context = ssl._create_unverified_context()

    next_url = BASE_URL
    attempts = 0
    output_lines = []

    while next_url:
        try:
            req = request.Request(next_url, headers={"Accept": "application/json"})
            res = request.urlopen(req, context=context)

            if res.status == 408:
                sleep(61)
                continue
            elif res.status == 204:
                break

            payload = json.loads(res.read().decode())
            next_url = payload["next"]
            attempts = 0

        except HTTPError as e:
            if e.code == 408:
                sleep(61)
                continue
            elif attempts < 3:
                attempts += 1
                sleep(61)
                continue
            else:
                sys.stderr.write(f"Failed for {ipr_id}. LAST URL: {next_url}\n")
                return

        for item in payload["results"]:
            entries = item.get("entries", None)

            if entries:
                entries_header = "-".join(
                    [entry["accession"] + "(" + ";".join(
                        [
                            ",".join(
                                [str(fragment["start"]) + "..." + str(fragment["end"]) 
                                 for fragment in locations["fragments"]]
                            ) for locations in entry["entry_protein_locations"]
                        ]
                    ) + ")" for entry in entries]
                )
                header = f">{item['metadata']['accession']}{HEADER_SEPARATOR}{entries_header}{HEADER_SEPARATOR}{item['metadata']['name']}"
            else:
                header = f">{item['metadata']['accession']}{HEADER_SEPARATOR}{item['metadata']['name']}"

            output_lines.append(header)

            seq = item["extra_fields"]["sequence"]
            for i in range(0, len(seq), LINE_LENGTH):
                output_lines.append(seq[i:i+LINE_LENGTH])

        if next_url:
            sleep(1)

    # Save to file
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, "w") as f:
        f.write("\n".join(output_lines) + "\n")
    print(f"[✓] Saved {ipr_id} to {filename}")
    print()

def main(ipr_ids, reviewed=True):
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        future_to_id = {
            executor.submit(fetch_sequences, ipr_id, reviewed): ipr_id
            for ipr_id in ipr_ids
        }
        for future in as_completed(future_to_id):
            ipr_id = future_to_id[future]
            try:
                future.result()
            except Exception as exc:
                sys.stderr.write(f"Error fetching {ipr_id}: {exc}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch InterPro FASTA sequences.")
    parser.add_argument("ipr_ids", nargs="+", help="List of IPR IDs to fetch")
    parser.add_argument("--unreviewed", action="store_true", help="Use unreviewed proteins instead of reviewed")

    args = parser.parse_args()
    main(args.ipr_ids, reviewed=not args.unreviewed)
