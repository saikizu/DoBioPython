# Author: Michael Do (saikizu)
# File purpose: Combine Domain and Sequence work
# Package: DoBioPython
# Date: 10-Jun-2024

import re


def Gap(letter: str, begin: int, end: int):
    """
    Create gap letters.\n
    For example: from 5 to 10 create gap with -\n
    Gap('-', 5, 10) -> "------"
    """
    gap = ""
    for i in range(begin - 1, end):
        gap += letter
    return gap


def Subseq(sequence: str, begin: int, end: int):
    '''
    Get the sub-sequence from 'begin' to 'end'.\n
    Example: trim seq A from 5 to 10.\n
    Subseq(seq, 5, 10) -> "ATCCTT"
    '''
    output = ""
    for i in range(begin - 1, end)
        output += sequence[i]
    return output


def CombineDomain(fileExtract: str, fileOutputExtract: str = None) -> dict:
    """ Combine all overlapping domain to one.

    Args:
        fileExtract (str): path of input file 3 first column must in this 
        format:\n
        "{gene} {from} {to}"
        fileOutputExtract (str, optional): save {path} text file as format:\n
        "{gene} {from} {to}". Defaults to None.

    Returns:
        dict: gene location dictionary from file
    """
    fileExtract = open(fileExtract, 'r').readlines()
    gene = {}
    for line in fileExtract:
        extract = re.match(r"(\S+)\s+(\d+)\s+(\d+)", line)
        try:
            loc = [int(extract.group(2)), int(extract.group(3))]
            if extract.group(1) in gene.keys():
                temp = gene[extract.group(1)]
                loc0 = ""
                loc1 = ""
                for i in range(len(temp)):
                    if loc[0] < temp[i][0]:
                        loc0 = f"<{i}"
                        break
                    elif loc[0] <= temp[i][1]:
                        loc0 = f"={i}"
                        break
                for i in range(len(temp)):
                    if loc[1] < temp[i][0]:
                        loc1 = f"<{i}"
                        break
                    elif loc[1] <= temp[i][1]:
                        loc1 = f"={i}"
                        break
                if loc0 + loc1 == "":
                    temp.append(loc)
                    gene[extract.group(1)] = temp
                    continue
                if loc0 != "":
                    start = int(loc0[1:])
                    loc0 = loc0[0]
                if loc1 != "":
                    end = int(loc1[1:])
                    loc1 = loc1[0]
                else:
                    end = len(temp)
                if start == end:
                    if loc0 + loc1 == "==":
                        continue
                    if loc0 + loc1 == "<<":
                        temp.insert(start, loc)
                        gene[extract.group(1)] = temp
                        continue
                    if loc0 + loc1 == "<=":
                        temp[start][0] = loc[0]
                        gene[extract.group(1)] = temp
                        continue
                else:
                    if loc0 + loc1 == "==":
                        temp[start][1] = temp[end][1]
                        del temp[end:start:-1]
                        gene[extract.group(1)] = temp
                        continue
                    if loc0 + loc1 == "<=":
                        temp[start][0] = loc[0]
                        temp[start][1] = temp[end][1]
                        del temp[end:start:-1]
                        gene[extract.group(1)] = temp
                        continue
                    if loc0 + loc1 == "=<":
                        temp[start][1] = loc[1]
                        del temp[end - 1:start:-1]
                        gene[extract.group(1)] = temp
                        continue
                    if loc0 + loc1 == "<<":
                        temp[start] = loc
                        del temp[end - 1:start:-1]
                        gene[extract.group(1)] = temp
                        continue
                    if loc0 + loc1 == "=":
                        temp[start][1] = loc[1]
                        if start < len(temp) - 1:
                            del temp[start + 1:]
                        gene[extract.group(1)] = temp
                        continue
                    if loc0 + loc1 == "<":
                        temp[start] = loc
                        if start < len(temp) - 1:
                            del temp[start + 1:]
                        gene[extract.group(1)] = temp
                        continue
            else:
                gene[extract.group(1)] = [loc]
        except TypeError:
            print(TypeError)
            continue
    if fileOutputExtract is not None:
        output = ""
        for key in gene.keys():
            for each in gene[key]:
                output += f"{key}\t{each[0]}\t{each[1]}\n"
        open(fileOutputExtract, "w").write(output)
    return gene


def ExtractDomainFromFasta(fileFasta: str, fileExtract: str,
                           gap_letter: str = None, fileOutputFasta: str = None,
                           fileOutputExtract: str = None) -> dict:
    """ Combine all overlapping domain to one. Then, extract all domain from
    fasta sequences.

    Args:
        fileFasta (str): path of FASTA file
        fileExtract (str): path of input file 3 first column must in this
        format:
        "{gene} {from} {to}"
        gap_letter (str, optional): The replacement for gaps. Defaults to None.
        fileOutputFasta (str, optional): save {path} output file as FASTA
        format. Defaults to None.
        fileOutputExtract (str, optional): save {path} output Domain as format:\n
        "{gene} {from} "{to}". Defaults to None.

    Returns:
        dict: gene sequences dictionary
    """
    gene = CombineDomain(fileExtract, fileOutputExtract)
    fileFasta = open(fileFasta, 'r').readlines()
    name = ""
    geneFasta = {}
    for line in fileFasta:
        if line.startswith(">"):
            name = line[1:].replace("\n", "")
            continue
        if bool(re.search(r"^\S", line)):
            try:
                if name in gene.keys():
                    subseq = []
                    for i in range(len(gene[name])+1):
                        try:
                            if gap_letter is not None:
                                if not subseq:
                                    subseq.append(
                                        Gap(gap_letter, 1,
                                            gene[name][i][0] - 1))
                                else:
                                    subseq.append(
                                        Gap(gap_letter,
                                            gene[name][i - 1][1] + 1,
                                            gene[name][i][0] - 1))
                            subseq.append(Subseq(line, gene[name][i][0],
                                                 gene[name][i][1]))
                        except IndexError:
                            subseq.append(Gap(gap_letter,
                                              gene[name][i - 1][1] + 1,
                                              len(line)))
                    geneFasta[name] = "".join(subseq)
                else:
                    print(f"{name} was not found")
            except IndexError:
                print(IndexError)
                continue
    if fileOutputFasta is not None:
        output = ""
        for key in geneFasta.keys():
            output += f">{key}\n{geneFasta[key]}\n"
        open(fileOutputFasta, 'w').write(output)
    return geneFasta
