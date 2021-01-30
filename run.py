import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import cosine_similarity
import os
import argparse


def get_feature(X, modelLen):
    def getVector(x, num, modelLen):
        vector = np.zeros(2 * modelLen + 1)
        if num <= modelLen:
            mean = x[0:num + 1].mean()
            vector[0:modelLen - num] = [mean] * (modelLen - num)
            vector[modelLen - num:len(vector)] = x[0:num + 1 + modelLen]
        elif num < (len(x) - modelLen):
            vector = x[num - modelLen:num + modelLen + 1]
        else:
            mean = x[i:].mean()
            vector[0:modelLen + (len(x) - num)] = x[num - modelLen:]
            vector[modelLen + (len(x) - num):] = [mean] * (modelLen + 1 + num - len(x))
        return vector

    s = np.zeros(shape=(len(X), 2 * modelLen + 1))
    for i in range(len(X)):
        vector = getVector(X[i], i, modelLen)
        s[i] = vector
    return s


def calLaplacianMatrix(adjacentMatrix):
    degreeMatrix = np.sum(adjacentMatrix, axis=1)
    laplacianMatrix = np.diag(degreeMatrix) - adjacentMatrix
    sqrtDegreeMatrix = np.diag(1.0 / (degreeMatrix ** 0.5 + 1e-5))
    return np.dot(np.dot(sqrtDegreeMatrix, laplacianMatrix), sqrtDegreeMatrix)


def get_domain(result):
    n = 0
    domain = []
    begin = 0
    end = 0
    try:
        while n < len(result):
            begin = n
            end = n
            temp = result[begin]
            while temp == result[n]:
                end = end + 1
                n = n + 1
            if (end - begin) >= 5:
                domain.append([begin, end - 1])
    except IndexError:
        if (end - 1 - begin) >= 5:
            domain.append([begin, end - 1])
    return domain


def write_domain(result_path, domain, chr_num):
    outfile = result_path + "domains.chr" + str(chr_num)
    f = open(outfile, 'w')
    domain_str = []
    for elem in domain:
        domain_str.append([str(elem[0]), str(elem[1])])
    for elem in domain_str:
        f.write(elem[0])
        f.write(",")
        f.write(elem[1])
        f.write("\n")
    f.close()


def write_domain_del(result_path, domain, chr_num):
    outfile = result_path + "domains_del.chr" + str(chr_num)
    f = open(outfile, 'w')
    domain_str = []
    for elem in domain:
        domain_str.append([str(elem[0]), str(elem[1])])
    for elem in domain_str:
        f.write(elem[0])
        f.write(",")
        f.write(elem[1])
        f.write("\n")
    f.close()


def write_dif(result_path, intra_list, inter_list, chr_num):
    outfile = result_path + "dif_diagnal.chr" + str(chr_num)
    f = open(outfile, 'w')
    difs = intra_list - inter_list
    f.write(str(difs.mean()) + "\n")
    for i in range(len(difs)):
        f.write('{:>6.3f}\t{:>6.3f}\t{:>6.3f}\t'.format(intra_list[i], inter_list[i], difs[i]))
        f.write("\n")
    f.close()


def intra(data, domains, diagnol=1):
    intra_list = np.zeros(len(domains))
    if diagnol:
        for i in range(len(domains)):
            temp = data[domains[i][0]:domains[i][1], domains[i][0]:domains[i][1]].mean()
            intra_list[i] = temp
    else:
        for i in range(len(domains)):
            tmp = []
            for j in range(domains[i][0], domains[i][1] + 1):
                for z in range(j + 1, domains[i][1] + 1):
                    tmp.append(data[j][z])

            intra_list[i] = np.mean(tmp)
    intra_average = intra_list.mean()
    return intra_list, intra_average


def inter(data, domains):
    inter_list = np.zeros(len(domains))
    for i in range(len(domains)):
        if i == 0:
            temp = data[domains[0][0]:domains[0][1], domains[1][0]:domains[1][1]].mean()

        elif i == (len(domains) - 1):
            temp = data[domains[i][0]:domains[i][1], domains[i - 1][0]:domains[i - 1][1]].mean()
        else:
            temp = (data[domains[i][0]:domains[i][1], domains[i - 1][0]:domains[i - 1][1]].mean() + data[domains[i][0]:
                    domains[i][1],domains[i + 1][0]:domains[i + 1][1]].mean()) / 2
        inter_list[i] = temp
    inter_average = inter_list.mean()
    return inter_list, inter_average


def SBTD(hicData, kmeans, outpath, chr_num):
    originData = hicData
    data_mean = originData.mean()
    data = np.log1p(hicData)
    features = get_feature(data, 50)
    sim = cosine_similarity(features) * 0.5 + 0.5
    Laplacian = calLaplacianMatrix(sim)
    w, v = np.linalg.eigh(Laplacian)
    H = v[:, 0: kmeans]
    sp_kmeans = KMeans(n_clusters=kmeans).fit(H)
    domains = get_domain(sp_kmeans.labels_)
    del_index = []
    for i in range(len(domains)):
        if originData[domains[i][0]:domains[i][1], domains[i][0]:domains[i][1]].mean() < data_mean:
            del_index.append(i)
    domains_del = [domains[i] for i in del_index]
    domains = np.delete(domains, del_index, 0)
    intra_list, intra_ave = intra(originData, domains)
    inter_list, inter_ave = inter(originData, domains)
    write_domain(outpath, domains, chr_num)
    write_domain_del(outpath, domains_del, chr_num)
    write_dif(outpath, intra_list, inter_list, chr_num)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input',type=str, help="Input file address")
    parser.add_argument("-k", '--cluster_Num', type=int, default=3,help="the number of clusters")
    parser.add_argument("-o",'--outpath',type=str, default="./output/",help="Output file address")
    parser.add_argument("-c","--chr_Num",type=int, default=1, help="Chromosome number")
    args = parser.parse_args()

    if not os.path.exists(args.outpath):
        os.mkdir(args.outpath)
    data = np.loadtxt(args.input)
    SBTD(data, args.cluster_Num, args.outpath, args.chr_Num)
