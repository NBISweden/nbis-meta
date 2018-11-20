import sys, pandas as pd, io

input_stream = io.TextIOWrapper(sys.stdin.buffer, encoding='utf-8', errors="ignore")

def readlines(fh):
    parsed = {}
    d = {}
    keys = ["AC","ID","DE","EC"]
    for i, line in enumerate(fh):
        line = line.rstrip()
        items = line.rsplit()
        try: key = items[0]
        except IndexError:
            continue
        if key in keys:
            if key == "ID" and i>0:
                ac = d["AC"]
                del(d["AC"])
                parsed[ac] = d
                d = {}
            d[key] = " ".join(items[1:])
    ac = d["AC"]
    del(d["AC"])
    parsed[ac] = d
    return parsed


def main():
    with input_stream as fh:
        parsed = readlines(fh)

    df = pd.DataFrame(parsed).T
    df = df[["ID", "DE", "EC"]]
    df.rename(columns={"ID": "TIGRFAM_ID","DE": "Description","EC": "Enzyme_ID"},inplace=True)
    df.index.name="Accession"
    df.to_csv(sys.stdout,sep="\t")

if __name__ == '__main__':
    main()