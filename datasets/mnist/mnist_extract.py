from array import array as pyarray
from scipy.io import loadmat
from sklearn.decomposition import PCA

import gzip
import hashlib
import logging
import numpy as np
import os
import os.path
import struct
import sys
import wget


TRAIN_IMAGES_URL = "http://yann.lecun.com/exdb/mnist/train-images-idx3-ubyte.gz"
TRAIN_LABELS_URL = "http://yann.lecun.com/exdb/mnist/train-labels-idx1-ubyte.gz"
TEST_IMAGES_URL  = "http://yann.lecun.com/exdb/mnist/t10k-images-idx3-ubyte.gz"
TEST_LABELS_URL  = "http://yann.lecun.com/exdb/mnist/t10k-labels-idx1-ubyte.gz"

TRAIN_IMAGES_SHA256 = "440fcabf73cc546fa21475e81ea370265605f56be210a4024d2ca8f203523609"
TRAIN_LABELS_SHA256 = "3552534a0a558bbed6aed32b30c495cca23d567ec52cac8be1a0730e8010255c"
TEST_IMAGES_SHA256 = "8d422c7b0a1c1c79245a5bcf07fe86e33eeafee792b84584aec276f5a2dbc4e6"
TEST_LABELS_SHA256 = "f7ae60f92e00ec6debd23a6088c31dbd2371eca3ffa0defaefb259924204aec6"

TRAIN_SAMPLE_INDICES_FNAME = "mnist_train_sample.tbl"
TEST_SAMPLE_INDICES_FNAME  = "mnist_test_sample.tbl"

FNAME_IMG = {
    'train': 'train-images-idx3-ubyte.gz',
    'test':  't10k-images-idx3-ubyte.gz'
}

FNAME_LBL = {
    'train': 'train-labels-idx1-ubyte.gz',
    'test':  't10k-labels-idx1-ubyte.gz'
}


def download_and_check(in_url, out_fname, sha256sum):
    logging.info("Downloading '{}'".format(in_url))
    wget.download(in_url, out_fname)

    valid = False
    with open(out_fname, "rb") as f:
        valid = (hashlib.sha256(f.read()).hexdigest() == sha256sum)

    return valid


def load_mnist(data="train", digits=np.arange(10)):
    fname_img = FNAME_IMG[data]
    fname_lbl = FNAME_LBL[data]

    with gzip.open(fname_lbl, 'rb') as flbl:
        magic_nr, size = struct.unpack(">II", flbl.read(8))
        lbl = pyarray("b", flbl.read())

    with gzip.open(fname_img, 'rb') as fimg:
        magic_nr, size, rows, cols = struct.unpack(">IIII", fimg.read(16))
        img = pyarray("B", fimg.read())

    ind = [k for k in range(size) if lbl[k] in digits]
    N = len(ind)

    images = np.zeros((N, rows*cols), dtype=np.uint8)
    labels = np.zeros((N, 1), dtype=np.int8)
    for i in range(len(ind)):
        m = ind[i]*rows*cols
        n = (ind[i]+1)*rows*cols
        images[i] = np.array(img[m:n])
        labels[i] = lbl[ind[i]]

    return images, labels


if __name__ == "__main__":
    logging.basicConfig(filename="mnist_extract.log",
                        format="%(levelname)s:%(message)s",
                        level=logging.INFO)

    # Get and check original data if needed
    urls = [TRAIN_IMAGES_URL, TRAIN_LABELS_URL,
            TEST_IMAGES_URL, TEST_LABELS_URL]
    fnames = [FNAME_IMG['train'], FNAME_LBL['train'],
              FNAME_IMG['test'], FNAME_LBL['test']]
    sha256sums = [TRAIN_IMAGES_SHA256, TRAIN_LABELS_SHA256,
                  TEST_IMAGES_SHA256, TEST_LABELS_SHA256]
    for url, fname, sha256sum in zip(urls, fnames, sha256sums):
        if not os.path.exists(fname):
            ok = download_and_check(url, fname, sha256sum)
            if not ok:
                logging.error("'{}' is corrupted; aborting".format(fname))
                exit(1)

    # We now have the original data
    logging.info("Loading MNIST training data")
    mnist_train = dict()
    mnist_train['train_X'], mnist_train['train_labels'] = load_mnist("train")
    train_size  = mnist_train['train_X'].shape[0]

    logging.info("Loading MNIST test data")
    mnist_test = dict()
    mnist_test['test_X'], mnist_test['test_labels'] = load_mnist("test")
    test_size  = mnist_test['test_X'].shape[0]

    should_load_samples = False
    if len(sys.argv) == 2 \
       or (not os.path.exists(TRAIN_SAMPLE_INDICES_FNAME)) \
       or (not os.path.exists(TEST_SAMPLE_INDICES_FNAME)):
        sample_size = int(sys.argv[1])

        if sample_size/2 > min(train_size, test_size):
            print("sample size is too large")
            should_load_samples = True
        else:
            logging.info("Generating {} samples".format(sample_size))
            train_sample_indices = np.randint(0, train_size, sample_size / 2)
            test_sample_indices  = np.randint(0, test_size, sample_size / 2)

            logging.info("Saving generated samples")
            np.savetxt("mnist_train_sample.tbl", train_sample_indices, fmt="%u")
            np.savetxt("mnist_test_sample.tbl", test_sample_indices, fmt="%u")
    else:
        should_load_samples = True

    if should_load_samples:
        logging.info("Loading samples")
        train_sample_indices = np.loadtxt(TRAIN_SAMPLE_INDICES_FNAME, dtype=int)
        test_sample_indices  = np.loadtxt(TEST_SAMPLE_INDICES_FNAME, dtype=int)
        sample_size = train_sample_indices.shape[0] \
                    + test_sample_indices.shape[0]

    logging.info("Extracting {} samples".format(sample_size))
    train_samples = mnist_train['train_X'][train_sample_indices, :]
    test_samples  = mnist_test['test_X'][test_sample_indices, :]
    mnist_sample  = np.concatenate((train_samples, test_samples))
    mnist_sample = PCA(n_components=512, whiten=True).fit_transform(mnist_sample)

    train_labels  = mnist_train['train_labels'][train_sample_indices]
    test_labels   = mnist_test['test_labels'][test_sample_indices]
    mnist_sample_labels = np.concatenate((train_labels, test_labels))

    logging.info("Saving extracted samples and their labels")
    sample_fname = "mnist_{}.tbl".format(sample_size)
    labels_fname = "mnist_{}.labels".format(sample_size)
    np.savetxt(sample_fname, mnist_sample, fmt="%f")
    np.savetxt(labels_fname, mnist_sample_labels, fmt="%u")
