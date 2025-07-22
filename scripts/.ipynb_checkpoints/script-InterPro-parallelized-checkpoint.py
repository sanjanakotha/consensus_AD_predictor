#!/usr/bin/env python3

import sys, json, ssl, os
from urllib import request
from urllib.error import HTTPError
from time import sleep
from concurrent.futures import ThreadPoolExecutor, as_completed

HEADER_SEPARATOR = "|"
LINE_LENGTH = 80
MAX_WORKERS = 5  # Adjust based on system/network capacity

def fetch_sequences(ipr_id):
    BASE_URL = f"https://www.ebi.ac.uk:443/interpro/api/protein/reviewed/entry/InterPro/{ipr_id}/?page_size=200&extra_fields=sequence"
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
    filename = f"../output/interpro/{ipr_id}.fasta"
    with open(filename, "w") as f:
        f.write("\n".join(output_lines) + "\n")
    print(f"[✓] Saved {ipr_id} to {filename}")

def main(ipr_ids):
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        future_to_id = {executor.submit(fetch_sequences, ipr_id): ipr_id for ipr_id in ipr_ids}
        for future in as_completed(future_to_id):
            ipr_id = future_to_id[future]
            try:
                future.result()
            except Exception as exc:
                sys.stderr.write(f"Error fetching {ipr_id}: {exc}\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write(f"Usage: {sys.argv[0]} <IPR_ID1> [<IPR_ID2> ...]\n")
        sys.exit(1)

    ipr_ids = sys.argv[1:]
    main(ipr_ids)
