from sklearn.decomposition import PCA
from sklearn.feature_extraction.text import TfidfVectorizer

import hashlib
import logging
import numpy as np
import os
import os.path
import sys
import tarfile
import wget


DATA_URL = "http://kdd.ics.uci.edu/databases/20newsgroups/20_newsgroups.tar.gz"
DATA_FILE = "20_newsgroups.tar.gz"
DATA_SHA256 = "b7bbf82b7831f7dbb1a09d9312f66fa78565c8de25526999b0d66f69d37e414"


def build_topic_corpus(corpus_file, n, topic):
    logging.info("Extracting corpus for topic '{}'".format(topic))
    topic_items = []
    names = corpus_file.getnames()
    for name in names:
        if topic in name:
            ti = corpus_file.getmember(name)
            if ti.isfile():
                topic_items.append(name)
    if len(topic_items) == 0:
        # Topic does not exist (no items fetched)
        raise ValueError(topic)

    topic_ids = []
    topic_corpus = []
    indices = np.arange(len(topic_items))
    np.random.shuffle(indices)
    indices = indices[:n]
    for i in indices:
        ti = corpus_file.getmember(topic_items[i])
        with corpus_file.extractfile(ti) as f:
            try:
                contents = str(f.read(), encoding="utf8")
            except ValueError as e:
                logging.warn("Encoding error in '{}': {}".format(ti.name, e))
                continue
            _, item_id = os.path.split(ti.name)
            topic_ids.append(item_id)
            topic_corpus.append(contents)

    return topic_ids, topic_corpus


def build_corpus(n, topics):
    """
    Builds a corpus with each topic, with N items each.
    Returns a list of document IDs and a corpus which is a dict where each topic
    is a key mapped to a list of document contents.
    """
    ids = []
    corpus = dict()
    with tarfile.open(DATA_FILE, "r:gz") as f:
        for topic in topics:
            topic_ids, topic_corpus = build_topic_corpus(f, n, topic)
            corpus[topic] = topic_corpus
            ids.extend(topic_ids)
    return ids, corpus


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("usage: {} STOP_WORDS N TOPIC [ TOPIC [ ... ] ]".format(sys.argv[0]))
        print("The program reads the file STOP_WORDS for stop words, extracts"
            + " and generates a BoW model from N random articles of each TOPIC")
        exit(1)

    logging.basicConfig(filename="newsgroups_extract.log",
                        format="%(levelname)s:%(message)s",
                        level=logging.INFO)

    if not os.path.exists(DATA_FILE):
        logging.info("Downloading data from '{}'".format(DATA_URL))
        wget.download(DATA_URL, DATA_FILE)
        with open(DATA_FILE, "rb") as f:
            if not hashlib.sha256(f.read()).hexdigest() != DATA_SHA256:
                logging.error("'{}' is corrupted; aborting".format(DATA_FILE))
                exit(1)

    # Read stop words list
    try:
        with open(sys.argv[1]) as stop_words_file:
            stop_words = stop_words_file.read().split()
    except Exception as e:
        logging.error("Could not read stop words: {}".format(e))
        exit(1)

    try:
        n = int(sys.argv[2])
        if (n < 2) or (n > 1000):
            raise ValueError("N must be between 2 and 1000")
    except ValueError as e:
        logging.error("Invalid argument: {}".format(e))
        exit(1)

    # Extract text corpus from tarball
    logging.info("Building corpus")
    topics = sys.argv[3:]
    try:
        ids, corpus = build_corpus(n, topics)
    except ValueError as e:
        logging.error("Invalid topic: {}".format(e))
        exit(1)

    corpus_text = []
    for topic_items in corpus.values():
        corpus_text.extend(topic_items)

    # Compute the TF-IDF matrix
    logging.info("Computing TF-IDF matrix")
    vectorizer = TfidfVectorizer(min_df=0.01, stop_words=stop_words)
    X = vectorizer.fit_transform(corpus_text)

    # Reduce data dimensionality using PCA
    logging.info("Computing PCA and reducing to 512 dimensions")
    X = PCA(n_components=512, whiten=True).fit_transform(X.toarray())

    # Save all extracted features and related data
    logging.info("Writing IDs file")
    ids_fname = "newsgroups-{}-{}.ids".format(n, len(topics))
    np.savetxt(ids_fname, ids, fmt="%s")

    logging.info("Writing table file")
    tbl_fname = "newsgroups-{}-{}.tbl".format(n, len(topics))
    np.savetxt(tbl_fname, X.todense(), fmt="%f")

    logging.info("Writing labels file")
    labels_fname = "newsgroups-{}-{}.labels".format(n, len(topics))
    counts = [len(topic_items) for topic_items in corpus.values()]
    np.savetxt(labels_fname, np.repeat(topics, counts), fmt="%s")
