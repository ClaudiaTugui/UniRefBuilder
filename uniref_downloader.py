import urllib
import urllib.request
import json
import os
import time
import pandas
from pathlib import Path
import shutil

cache_path = 'uniref_cache'
if not os.path.exists(cache_path):
    print(f"./{cache_path}, does not exist, creating...")
    os.makedirs(cache_path)
    print(f"./{cache_path} created")
else:
    print(f"./{cache_path} already exists")


def retry_until_successful(work):
    success = False
    while not success:
        try:
            return work()
        except KeyboardInterrupt:
                raise KeyboardInterrupt
        except Exception as e:
                print(e)
                print("Something went wrong with the request, retrying")
        time.sleep(5)

def call_uniprot(url):
    response = urllib.request.urlopen(url, timeout=10)
    response_text = response.read().decode('UTF-8')
    results = json.loads(response_text)
    total_records = int(response.headers["x-total-records"])
    cursor =  response.headers["link"].split(" ")[0].replace("<", "").replace(">", "")[:-1] if response.headers["link"] else None
    return [results["results"], total_records, cursor]

def download_uniref_data_for_organism(organism, identity):
    print(f"Downloading data for {organism}...")

    [i, cursor] = [0, f"https://rest.uniprot.org/beta/uniref/search?query=({organism})%20AND%20(identity:{identity})&size=500"]
    cursor_filename = f"{cache_path}/{organism}.cursor"
    try:
        with open(cursor_filename, "r") as cursor_file:
            cursor_file_contents = cursor_file.read()
            [i, cursor] = cursor_file_contents.split("\t")
            i = int(i)
            print(f"cursor has been cached from last run, continuing from index {i}")
    except FileNotFoundError :
        print(f"cursor has not been cached, starting from index {i}")
        pass

    while cursor and cursor != "None":
        fasta_file_conents = ""

        print(cursor)
        [results, total_records, cursor] = retry_until_successful(lambda: call_uniprot(cursor))

        print(f"downloaded {round(((100/(total_records))*(len(results)+500*(i))), 3)}% of {organism} records ({len(results)+500*(i)} of {total_records})")

        for result in results:
            fasta_entry = f">{result['id']} {result['representativeMember']['proteinName']} n={result['organismCount']} Tax={result['commonTaxon']['scientificName']} TaxID={result['commonTaxon']['taxonId']} RepID={result['representativeMember']['memberId']}"

            fasta_entry += "\n"
            fasta_entry += result["representativeMember"]["sequence"]["value"]

            fasta_file_conents += fasta_entry
            fasta_file_conents += "\n\n"

        fasta_filename = f"{cache_path}/{organism}_{i}.fasta"
        with open(fasta_filename, "w") as fasta_file:
            fasta_file.write(fasta_file_conents)
        with open(cursor_filename, "w") as cursor_file:
            cursor_file.write(f"{i}\t{cursor}")

        i += 1

def FASTA_merger(clean_up=False):
    size_dir = 0
    database_file = "database.fasta.final"
    if os.path.exists(database_file):
        os.remove(database_file)
    files = Path(cache_path).glob("*.fasta")
    for file in files:
        with open(file) as f:
            FASTA = f.read()
            final = open(database_file, "a")
            final.write("\n")
            final.write(FASTA)
            final.close()
            size_dir += os.path.getsize(file)
    
    # print("Size of all files : "+str(size_dir) +" and final fasta file : "+ \
    #               str(os.path.getsize(database_file)))
    if clean_up == True:
        if abs(size_dir - os.path.getsize(database_file)) < 500:
            shutil.rmtree(cache_path)
            print("Deleted storage directory")
        else:
            print("Size of final file and seperate entries does not match, did not delete")

csv_filename = input("Please provide .csv file with genus names : ")
print("Identity options : ")
print("\tUNIREF 50 = 0.5")
print("\tUNIREF 90 = 0.9")
print("\tUNIREF 100 = 1.0")
identity = None
supported_identities =  ["0.5", "0.9", "1.0"]
while identity not in supported_identities:
    identity = input("Please input your identity option : ")
    if identity not in supported_identities:
        print("Unsupported option, please try again.")

csv_contents = pandas.read_csv(csv_filename, header=None)
organisms = csv_contents[0]
for organism in organisms:
    download_uniref_data_for_organism(organism, identity)
    FASTA_merger(clean_up=False)
