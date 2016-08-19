from scipy.io import loadmat
from scipy.misc import imsave
from sklearn.decomposition import PCA

import hashlib
import logging
import numpy as np
import os
import os.path
import sklearn.decomposition
import subprocess
import wget


# Original data
DATA_URL      = "http://isomap.stanford.edu/face_data.mat.Z"
SHA256_DIGEST = "9c5bc75f204071bbd340aa3ff584757ec784b0630206e526d4cd3809f2650a8a"

# Local name
DATA_FNAME = "face_data.mat"

# Output files/directories
IMG_DIR      = "images"
IMG_FNAME    = "face_raw.tbl"
LIGHTS_FNAME = "face_lights.tbl"
POSES_FNAME  = "face_poses.tbl"
PCA_FNAME    = "faces.tbl"


if __name__ == "__main__":
    logging.basicConfig(filename="faces_extract.log",
                        format="%(levelname)s:%(message)s",
                        level=logging.INFO)

    # Get original data
    if not os.path.exists(DATA_FNAME):
        if not os.path.exists("{}.Z".format(DATA_FNAME)):
            logging.info("Downloading faces data from '{}'".format(DATA_URL))
            wget.download(DATA_URL, "{}.Z".format(DATA_FNAME))

        logging.info("Checking SHA-1 digest")
        with open("{}.Z".format(DATA_FNAME), "rb") as f:
            if hashlib.sha256(f.read()).hexdigest() != SHA256_DIGEST:
                logging.error("File seems corrupted; aborting")
                exit(1)

        logging.info("Uncompressing data into '{}'".format(DATA_FNAME))
        subprocess.call(["uncompress", "{}.Z".format(DATA_FNAME)])

    # We have the original data; proceed
    logging.info("Loading faces data")
    faces = loadmat(DATA_FNAME)

    face_images = faces["images"]
    logging.info("Writing image table data to {}".format(IMG_FNAME))
    np.savetxt(IMG_FNAME, face_images.T, fmt="%f")

    if not os.path.exists(IMG_DIR):
        logging.info("Creating directory {}".format(IMG_DIR))
        os.makedirs(IMG_DIR, 0o755)
    elif not os.path.isdir(IMG_DIR):
        logging.error("File {} exists; aborting".format(IMG_DIR))
        exit(1)

    logging.info("Writing image files to {}".format(IMG_DIR))
    for i in range(face_images.shape[1]):
        image = face_images[:, i]
        image = image.reshape(64, 64).T
        path = os.path.join(IMG_DIR, "{}.png".format(i))
        imsave(path, image)

    logging.info("Writing lights data to {}".format(LIGHTS_FNAME))
    np.savetxt(LIGHTS_FNAME, faces["lights"].T, fmt="%f")

    logging.info("Writing poses data to {}".format(POSES_FNAME))
    np.savetxt(POSES_FNAME, faces["poses"].T, fmt="%f")

    logging.info("Writing PCA-whitened data to {}".format(PCA_FNAME))
    X = faces["images"].T
    X = PCA(n_components=256, whiten=True).fit_transform(X)
    np.savetxt(PCA_FNAME, X, fmt="%f")
