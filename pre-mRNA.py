__author__ = 'Guoliang Lin'
def trim(mdata):
    pre=str(mdata)
    pre=pre.replace(' ','')
    pre=pre.replace('\'','')
    pre=pre.replace('\\n','')
    pre=pre.replace('[','')
    pre=pre.replace(']','')
    pre=pre.replace(';','')
    pre=pre.replace('"','')
    pre=pre.replace(',','\t')
    return pre

def massgenerat(chrom):
    with open(listname[0][-1], 'w') as outfile:
        with open(listname[1][-1], 'w') as outfiler:
            with open(listname[2][-1], 'w') as outfiled:
                for each_line in infile:
                    list = each_line.split('\t')
                    tmp=list.pop().split(' ')
                    list.append(tmp)
                    if list[2]=='transcript':
                        if list[0] == chrom:
                            if list[6] == "+":
                                outfile.write(trim(list) + '\n')
                            elif list[6] == '-':
                                outfiler.write(trim(list) + '\n')
                            else:
                                outfiled.write(trim(list) + '\n')
    infile.seek(0)

def addname(x):
    listname[0].append("chrom" + x + '.cls')
    listname[1].append("chrom_reverse" + x + '.cls')
    listname[2].append("chrom_dot" + x + '.cls')

list1 = []
listname = [[], [], []]
MaxMin = [[0], [10 ** 30]]
isoforms = []
tmp = ''
transcript = []
gene = []
with open("transcript.gtf") as infile:
    for number in range(1, 20):
        strnmb = str(number)
        addname(strnmb)
        massgenerat(strnmb)
    infile.seek(0)
    addname('X')
    massgenerat('X')
    addname('Y')
    massgenerat('Y')
    addname('MT')
    massgenerat('MT')
with open("gene_ID.data", 'w') as gene_ID:
    with open("transcript_ID", 'w') as transcript_ID:
        with open("findsplice.data", 'w') as findsplice:
            for filenamelist in listname:
                for filename in filenamelist:
                    with open(filename) as file:
                        isoforms.append(file.readline())
                        while isoforms[-1]:
                            print(isoforms[-1])
                            isofromsList = isoforms[-1].split('\t')
                            MaxMin[0].append(int(isofromsList[4]))
                            MaxMin[1].append(int(isofromsList[3]))
                            if ((MaxMin[0][0] >= MaxMin[0][-1]) & (MaxMin[1][0] <= MaxMin[0][-1])) | (
                                        (MaxMin[0][0] >= MaxMin[1][-1]) & (MaxMin[1][0] <= MaxMin[0][-1])) | (
                                        (MaxMin[0][0] < MaxMin[0][-1]) & (MaxMin[1][0] > MaxMin[1][-1])):
                                MaxMin[0][0] = max(MaxMin[0][0], MaxMin[0][-1])
                                MaxMin[1][0] = min(MaxMin[1][0], MaxMin[1][-1])
                            elif len(isoforms) > 2:
                                tmp = isoforms.pop()
                                findsplice.write(str(len(isoforms)) + '\n')
                                for item in isoforms:
                                    list1 = item.split('\t')
                                    if not ('CUFF' in list1[9]):
                                        gene_ID.write(list1[9] + '\n')
                                    if not ('CUFF' in list1[11]):
                                        transcript_ID.write(list1[11] + '\n')
                                    findsplice.write(item)
                                findsplice.write('\n')
                                MaxMin[0] = [MaxMin[0][-1]]
                                MaxMin[-1] = [MaxMin[1][-1]]
                                isoforms = [tmp]
                            else:
                                tmp = isoforms.pop()
                                MaxMin[0] = [MaxMin[0][-1]]
                                MaxMin[-1] = [MaxMin[1][-1]]
                                isoforms = [tmp]
                            isoforms.append(file.readline())
                        isoforms = []
                        MaxMin = [[0], [10 ** 30]]