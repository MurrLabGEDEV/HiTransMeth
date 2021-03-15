#!/usr/bin/env python3
# -*- coding: utf8 -*-

"""
:Author: Victor Ythier
:Date: 06/09/2019

Python module containing basic and common functions for the HiTransMeth Pipeline
    - setSemaphores
    - getSemaphores
    - check_filePath
    - check_fileSize
    - check_fileRefs
    - copy_wBackup

"""
import os, shutil

def setSemaphore(filename, value):
    """
    Set a semaphore file to a given value. The value has to be DONE, RUNNING, ERROR else an error is raised.

    :param filename: Path of the semaphore file.
    :param value: Value to be written (DONE, RUNNING, ERROR)
    :return: None if everything went fine
    """
    if value not in ['DONE', 'RUNNING', 'ERROR']:
        raise RuntimeError('A semaphore value can only be DONE, RUNNING or ERROR')

    with open(filename, 'w') as f:
        f.write(value)


def getSemaphore(filename):
    """
    Reads the semaphore file and return the value

    :param filename: Path of the semaphore file.
    :return: The semaphore value (DONE, RUNNING, ERROR) or None if the semaphore file does not exists.
    """
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            return f.readline().strip()
    else:
        return None


def check_filePath(files):
    """
    Checks that all files exist.

    :param files: list of files to check
    :return:
    """
    missing_files = []
    for fn in files:
        if not os.path.exists(fn):
            missing_files.append(fn)
    if len(missing_files) > 0:
        raise RuntimeError('Some files are missing (see below). Check their path!\n' +
                           '\n'.join(missing_files))
    return True


def check_fileSize(files, minsize):
    """
    Check if there exists a file whose size is below the threshold.

    :param files: list of files
    :param minsize: minimum size of file in kb
    :return: true if all files are bigger than 'minsize'. Raises a RuntimeError otherwise
    """
    small_files = []
    for fn in files:
        if not os.path.isdir(fn) and os.path.getsize(fn) < minsize*1024:
            small_files.append(fn)
    if len(small_files) > 0:
        raise RuntimeError('Some files are smaller than ' + str(minsize) + 'kb (see below).\n' + '\n'.join(small_files))
    return True


def check_fileRefs(ref, motif):
    """
    Check if there exists a motif reference folder for every motif in motif table give in the CONFIG_FILE.

    :param ref: Dictionary of reference {'motif_name': "path/to/ref/folder"}
    :param motif: list of motif
    :return: true if all motif folder exist in the reference folder
    """
    missing_motif = []
    for m in motif:
        if not os.path.isdir(ref[m]):
            missing_motif.append(ref)
    if set(ref.keys()) != set(motif)
        raise RuntimeError('Some references motif are missing (see below).\n' + '\n'.join(missing_motif))
    return True


def copy_wBackup(orig, dest):
    """
    Copy a file into a directory. If the file already exists in the directory, then move it into a *.bkpN
    filename where N is an increment.

    :param orig: Absolute path of the file to copy
    :param dest: Absolute path of the directory where the file is to be copied
    :return: None
    """
    if os.path.isdir(dest):
        raise RuntimeError(dest + " is a directory. A file was expected")

    if os.path.exists(dest):
        i = 1
        bkp = dest + '.bkp' + str(i)
        while os.path.exists(bkp):
            i += 1
            bkp = dest + '.bkp' + str(i)
        shutil.move(dest, bkp)
    shutil.copyfile(orig, dest)