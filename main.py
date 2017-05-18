

import argparse
import math
import os
import random
import subprocess

from collections import defaultdict
from datetime import datetime

class NadphPheno:
    def __init__(self, filename):
        self.rows   = []
        self.values = []

        with open(filename, 'r') as infile:
            for line in infile:
                cols = line.split()

                if len(cols) != 3:
                    print("skipping rows \"", line, "\"")
                    continue

                self.rows.append(cols[:-1])
                self.values.append(cols[2])

    def generateShuffledFile(self, resultName):

        random.shuffle(self.values)

        with open(resultName, 'w') as outfile:
            for i in range(len(self.rows)):
                row = self.rows[i]
                value = self.values[i]

                outfile.write('%s\t%s\t%s\n' % (row[0], row[1], value))


class MartViewType:
    def __init__(self, cols):
        self.geneName = cols[0]
        self.hashId = int(cols[0][4:])
        self.min = int(cols[1])
        self.max = int(cols[2])
        self.chromosome = cols[3]

class MartViewBiTree:
    class Node:
        def __init__(self):
            self.min = 0
            self.max = 0
            self.values = []

    def init(self, data):

        # increases max and decreases min by this value to allow more SNPS
        # to be matched with these genes
        adjustRange = 5000

        self.min = data[0].min
        self.max = data[0].max

        for gene in data:
            # dont allow gene.min to be negative after adjustment
            gene.min = max(0, gene.min - adjustRange)
            gene.max += adjustRange

            self.min = min(gene.min, self.min)
            self.max = max(gene.max, self.max)

        self.diff = self.max - self.min
        self.nodes = []
        self.num = len(data)

        sortedData = sorted(data, key=lambda x: x.min)

        inRange = []
        inRangeLowest = 0

        for elem in sortedData:
            if len(inRange) == 0:
                inRange.append(elem)
                inRangeLowest = elem.min
                continue

            if inRange[-1].max < elem.min:
                node = self.Node()

                node.min = inRangeLowest
                node.max = inRange[-1].max
                node.values = inRange

                self.nodes.append(node)

                inRange = [elem]
                inRangeLowest = elem.min
                continue

            if inRange[0].max < elem.min and len(inRange) > 3:
                node = self.Node()

                node.min = inRangeLowest
                node.max = elem.min - 1
                node.values = inRange

                self.nodes.append(node)

                for i, e in enumerate(inRange):
                    if e.max >= elem.min:
                        inRange = inRange[i:]
                        break

                inRangeLowest = elem.min

                inRange.append(elem)
                inRange = sorted(inRange, key=lambda x: x.max)
                continue

            inRange.append(elem)
            inRange = sorted(inRange, key=lambda x: x.max)

        if len(inRange) > 0:
            node = self.Node()

            node.min = inRangeLowest
            node.max = inRange[-1].max
            node.values = inRange

            self.nodes.append(node)

        self.numNodes = len(self.nodes)

    def find(self, value):
        if value < self.min or value > self.max:
            return []

        first = 0
        last = len(self.nodes)

        while first <= last:
            middle = (first + last) // 2

            node = self.nodes[middle]

            if value < node.min:
                last = middle - 1
            elif value > node.max:
                first = middle + 1
            else:
                return node.values

        return []

class MartView:

    def __init__(self, filename):
        self.trees = defaultdict(MartViewBiTree)

        values = defaultdict(list)

        with open(filename, 'r') as infile:
            counter = 0
            for line in infile:
                cols = [x.strip() for x in line.split(',')]

                try:
                    value = MartViewType(cols)
                    values[value.chromosome].append(value)
                except ValueError:
                    if counter != 0:
                        print("MartView invalid integer, ", cols)
                        exit(-1)

                    counter += 1
                    continue


        for key, data in values.items():
            self.trees[key].init(data)

    def findGenes(self, chromosome, value):
        for v in self.trees[chromosome].find(value):
            if value >= v.min and value <= v.max:
                yield v

# Main ######################################################################################################################

def entry():

    parser = argparse.ArgumentParser()
    parser.add_argument('-start', type=int, default=0)
    parser.add_argument('-end', type=int, default=10000)
    args = parser.parse_args()

    nadphPheno = NadphPheno('nadph_pheno_men.txt')
    martView = MartView('mart_export.txt')

    for i in range(args.start, args.end):

        resultDirectory = 'worker/results%d/' % i
        if not os.path.exists(resultDirectory):
            os.makedirs(resultDirectory)

        nadphPheno.generateShuffledFile(resultDirectory + 'nadph_pheno.txt')

        timePlinkStart = datetime.now()

        p = subprocess.Popen(
            [
                './plink',
                '--bfile', 'dgrp_filtered',
                '--covar', 'nadph_covariate.txt',
                '--covar-number', '1-11',
                '--linear',
                '--maf', '0.02',
                '--mpheno', '1',
                '--out',   resultDirectory + 'result',
                '--pheno', resultDirectory + 'nadph_pheno.txt'
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)

        out, err = p.communicate()

        if p.returncode != 0:
            print(out.decode('ascii', errors='surrogateescape'))
            print(err.decode('ascii', errors='surrogateescape'))
            break

        timePlinkEnd = datetime.now()

        # Step 2 ############################################################################################################

        resultFilename = resultDirectory + 'result.assoc.linear'

        geneTStatDict = defaultdict(float)
        naCounter = 0

        with open(resultFilename, 'r') as infile:
            for line in infile:
                cols = line.split()

                if cols[4] != 'ADD':
                    continue

                col0value = 'ERROR'

                if cols[0] == '1':
                    col0value = '2L'
                elif cols[0] == '2':
                    col0value = '2R'
                elif cols[0] == '3':
                    col0value = '3L'
                elif cols[0] == '4':
                    col0value = '3R'
                elif cols[0] == '5':
                    col0value = 'X'
                elif cols[0] == '6':
                    col0value = '4'
                else:
                    continue

                cols = [col0value, cols[1], cols[2], cols[8]]

                for gene in martView.findGenes(col0value, int(cols[2])):
                    if cols[-1] == 'NA': # TODO as if this is accurate, skip NA values?
                        naCounter += 1
                        continue
                    geneTStatDict[gene.hashId] += math.log10(float(cols[-1]))

        timeStep2End = datetime.now()

        os.remove(resultFilename)

        with open(resultDirectory + 'tstat.txt', 'w') as outfile:
            for key, tvalue in geneTStatDict.items():
                outfile.write('%s\t%f\n' % (key, -2.0 * tvalue))

        with open(resultDirectory + 'done.txt', 'w') as outfile:
            outfile.write('plink: ' + str(timePlinkEnd - timePlinkStart) + '\n')
            outfile.write('step2/3: ' + str(timeStep2End - timePlinkEnd) + '\n')
            outfile.write('NA Count: ' + str(naCounter) + '\n')

# Script Entry ##############################################################################################################

entry()
