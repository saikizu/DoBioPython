import re

# Input information before run
fileFasta = "Australasica_NLR.fasta"
fileExtract = "Aus_TIRLRR.txt"
# Save file
saveExtract = True
fileOutputExtract = "OutputExtract.txt"
saveFasta = True
fileOutputFasta = "OutputFasta.txt"

#Main script
fileExtract = open(fileExtract, 'r').readlines()
gene = {}
for line in fileExtract:
    extract = re.match(r"([^\t]+)\t([^\t]+)\t([^\t]+)", line)
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
if saveExtract:
    output = ""
    for key in gene.keys():
        for each in gene[key]:
            output += f"{key}\t{each[0]}\t{each[1]}\n"
    open(fileOutputExtract, "w").write(output)
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
                for each in gene[name]:
                    subseq.append(line[each[0]+1:each[1]+2])
                geneFasta[name] = "".join(subseq)
            else:
                print(f"{name} was not found")
        except IndexError:
            print(IndexError)
            continue
if saveFasta:
    output = ""
    for key in geneFasta.keys():
        output += f">{key}\n{geneFasta[key]}\n"
    open(fileOutputFasta, 'w').write(output)
