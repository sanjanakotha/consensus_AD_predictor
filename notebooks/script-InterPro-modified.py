#!/usr/bin/env python3

import sys, errno, re, json, ssl
from urllib import request
from urllib.error import HTTPError
from time import sleep

HEADER_SEPARATOR = "|"
LINE_LENGTH = 80

def output_list(ipr_id):
    BASE_URL = f"https://www.ebi.ac.uk:443/interpro/api/protein/reviewed/entry/InterPro/{ipr_id}/?page_size=200&extra_fields=sequence"

    # Disable SSL verification to avoid config issues
    context = ssl._create_unverified_context()

    next_url = BASE_URL
    last_page = False
    attempts = 0

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
            if not next_url:
                last_page = True

        except HTTPError as e:
            if e.code == 408:
                sleep(61)
                continue
            elif attempts < 3:
                attempts += 1
                sleep(61)
                continue
            else:
                sys.stderr.write("LAST URL: " + next_url + "\n")
                raise e

        for item in payload["results"]:
            entries = item.get("entries", None)

            if entries is not None:
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
                sys.stdout.write(">" + item["metadata"]["accession"] + HEADER_SEPARATOR
                                  + entries_header + HEADER_SEPARATOR
                                  + item["metadata"]["name"] + "\n")
            else:
                sys.stdout.write(">" + item["metadata"]["accession"] + HEADER_SEPARATOR 
                                 + item["metadata"]["name"] + "\n")

            seq = item["extra_fields"]["sequence"]
            fastaSeqFragments = [seq[i:i+LINE_LENGTH] for i in range(0, len(seq), LINE_LENGTH)]
            for fragment in fastaSeqFragments:
                sys.stdout.write(fragment + "\n")

        if next_url:
            sleep(1)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: {} <IPR_ID>\n".format(sys.argv[0]))
        sys.exit(1)

    ipr_id = sys.argv[1]
    output_list(ipr_id)
