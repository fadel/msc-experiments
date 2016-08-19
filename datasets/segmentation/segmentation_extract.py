import hashlib
import logging
import pandas as pd
import os
import os.path
import wget


DATA_URL = "https://archive.ics.uci.edu/ml/machine-learning-databases/image/segmentation.test"
DATA_SHA256 = "2e9e966479d54c6aaec309059376dd9c89c1b46bf3a23aceeefb36d20d93a189"
DATA_FILE = "segmentation.test"


if __name__ == "__main__":
    logging.basicConfig(filename="segmentation_extract.log",
                        format="%(levelname)s:%(message)s",
                        level=logging.INFO)

    if not os.path.exists(DATA_FILE):
        logging.info("Downloading '{}'".format(DATA_URL))
        wget.download(DATA_URL, DATA_FILE)
        with open(DATA_FILE, "rb") as f:
            if hashlib.sha256(f.read()).hexdigest() != DATA_SHA256:
                logging.error("{} is corrupted; aborting".format(DATA_FILE))


    df = pd.read_table(DATA_FILE, header=None, skiprows=4, delimiter=",")

    # First column contains class names, which we convert to numbers using the
    # 'class_labels' dict
    classes = set(df[0])
    numbers = [i for i in range(len(classes))]
    class_labels = dict(zip(classes, numbers))

    data = df.drop([0, 3], axis=1)
    data.to_csv("segmentation.tbl", sep=" ", index=False, header=False)

    labels = df[0].apply(lambda x: class_labels[x])
    labels.to_csv("segmentation.labels", sep=" ", index=False, header=False)
