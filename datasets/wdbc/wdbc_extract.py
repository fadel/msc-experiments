import hashlib
import logging
import pandas as pd
import os
import os.path
import wget


DATA_URL = "http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
DATA_SHA256 = "d606af411f3e5be8a317a5a8b652b425aaf0ff38ca683d5327ffff94c3695f4a"
DATA_FILE = "wdbc.data"


if __name__ == "__main__":
    logging.basicConfig(filename="wdbc_extract.log",
                        format="%(levelname)s:%(message)s",
                        level=logging.INFO)

    if not os.path.exists(DATA_FILE):
        logging.info("Downloading '{}".format(DATA_URL))
        wget.download(DATA_URL, DATA_FILE)
        with open(DATA_FILE, "rb") as f:
            if hashlib.sha256(f.read()).hexdigest() != DATA_SHA256:
                logging.error("'{}' is corrupted; aborting".format(DATA_FILE))
                exit(1)

    data = pd.read_table(DATA_FILE, header=None, delimiter=",")
    wdbc_ids = data[0]
    wdbc_labels = data[1]
    wdbc = data.drop([0, 1], axis=1)

    wdbc.to_csv("wdbc.tbl", sep=" ", index=False, header=False)
    wdbc_labels.to_csv("wdbc.labels", sep=" ", index=False, header=False)
    wdbc_ids.to_csv("wdbc.ids", sep=" ", index=False, header=False)
