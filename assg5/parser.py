##-- encoding=utf-8
##-- written by twins in NCI
##-- to run train file to get PCFG
import sys
import json
from collections import defaultdict
import types

#将单词翻译为相应的词
class PCFGParser():

    def __init__(self, modelname):
        self.modelname = modelname
        self.model = self.readMode(self.modelname)
        self.unaries = self.getUnaryRule()
        self.binaries = self.getBinarys()

    def writeOutResult(self, s, output=''):
        result, value = self.parse(s)
        stataics = self.writeBridge(s)
        f = file(output, 'w')
        f.write(result + '\n')
        print result
        f.write('%s\n' %value)
        print value
        for i in stataics:
            f.write('%s\n' %i)
            print i
        f.close
    #根据PCFG反推出相应的结果树
    def parse(self, s):
        temp = self.CKY(s.split(' '))
        result = json.dumps(temp[0]).replace('[','(').replace(']',')').replace(',','').replace('"','')
        return result, temp[1]
           
    #获得所有的非终结符
    def getNonterminal(self):
        models = self.model
        li = []
        for i in models:
            li.append(i[0])
        return li
    
    #查找X->w对应的概率
    def query_w(self, X, w):
        s = self.unaries
        if tuple([X, w]) in s.keys():
            return s.get(tuple([X, w]))
        else:
            return 0

    #查找A->BC的概率
    def query_x_y_z(self, x, y, z):
        s = self.binaries
        if tuple([x, y, z]) in s.keys():
            return s.get(tuple([x, y, z]))
        else:
            return 0
    
    #获取非终结符到终结符的元组映射
    def getUnaryRule(self):
        models = self.model
        rules = {}
        for i in models:
            temp = i[1].split()
            if len(temp) == 1:
                rules[tuple([i[0], i[1]])] = i[2]
        return rules

    #获取非终结符到终结符的映射
    def getBinarys(self):
        models = self.model
        names = {}
        for i in models:
            temp = i[1].split()
            if len(temp) == 2:
                names[tuple([i[0], temp[0], temp[1]])] = i[2]
        return names

    #读取训练好的模型文件
    def readMode(self, filename=''):
        f = file(filename, 'r')
        li = []
        while True:
            line = f.readline()
            if len(line) == 0:
                break
            lines = line.replace('\n', '').split('#')
            li.append((lines[0], lines[1], float(lines[2])))
        f.close() # close the file
        return li

    def CKY(self, x):
        '''
            CKY算法和viterbi算法的综合求最大可能parser。
            1. 训练相应的PCFG，然后找到可能的分析树
            2. 计算W(i, j)概率最大的概率
        '''
        n = len(x) 
        pi = defaultdict(float) 
        bp = {} 
        N = self.getNonterminal() 
        for i in xrange(n):
            w = x[i] 
            for X in N:
               pi[i, i, X] = self.query_w(X, w) #chart(j,j)
        for l in xrange(1, n): 
            for i in xrange(n-l):
                j = i + l
                for X in N:
                    max_score = 0
                    args = None
                    for R in self.binaries.keys(): 
                        if R[0] == X: 
                            Y, Z = R[1:]
                            for s in xrange(i, j):
                                if pi[i, s, Y] and pi[s + 1, j, Z]: # calculate score if both pi entries have non-zero score
                                    score = self.query_x_y_z(X, Y, Z) * pi[i, s, Y] * pi[s + 1, j, Z]
                                    if max_score < score:
                                        max_score = score
                                        args = Y, Z, s
                    if max_score: # update DP table and back pointers
                        pi[i, j, X] = max_score
                        bp[i, j, X] = args
        if pi[0, n-1, 'S']:
            return self.recover_tree(x, bp, 0, n-1, 'S'), pi[0, n-1, 'S']
        else: # if the tree does not have the start symbol 'S' as the root
            max_score = 0
            args = None
            for X in N:
                if max_score < pi[0, n-1, X]:
                    max_score = pi[0, n-1, X]
                    args = 0, n-1, X
            return self.recover_tree(x, bp, *args)

    #将内向距离和外向距离组合成字串
    def writeBridge(self, y):
        output = []
        inside, outside = self.insidePro(y), self.outsidePro(y, self.insidePro(y))
        for i in inside.keys():
            for j in outside.keys():
                #中间节点值相同
                if i[1] == j[1]:
                    if i[0] == j[0] and i[2] == j[2] and (inside.get(i) != 0 or outside.get(j) != 0):
                        s = '%s # %s # %s # %s # %s' %(j[1], j[0], j[2], inside.get(i), outside.get(j))
                        output.append(s) 
        return output


    #计算外向最大概率，基于CYK算法
    def outsidePro(self, y, inside=defaultdict(float)):
        x = y.split()
        n = len(x) 
        outside = defaultdict(float)
        N = self.getNonterminal() 
        outside[1, 'S', n] = 1
        index = [n - i for i in range(n-1)]
        for l in index:
            for s in xrange(1, n-l+2):
                for t in xrange(1, l):
                    p1 = s
                    e1 = s + l - 1
                    q1 = s + t - 1
                    p2 = s + t
                    e2 = s
                    q2 = s + l
                    for A in N:
                        for R in self.binaries.keys():
                            if R[0] == A:
                                B, C = R[1:]
                                if outside[p1, A, e1] and inside[q1+1, C, e1]:
                                    outside[p1, B, q1] = outside[p1, B, q1] + outside[p1, A, e1] * inside[q1+1, C, e1] * self.query_x_y_z(A, B, C)
                                if outside[e2, A, q2] and inside[e2, B, p2-1]:
                                    outside[p2, C, q2] = outside[p2, C, q2] + outside[e2, A, q2] * inside[e2, B, p2-1] * self.query_x_y_z(A, B, C)
        return outside

    #计算内向概率值,基于CYK算法
    def insidePro(self, y):
        x = y.split()
        n = len(x) 
        inside = defaultdict(float)
        N = self.getNonterminal() 
        for k in xrange(1, n+1):
            w = x[k-1] 
            for A in N:
               inside[k, A, k] = self.query_w(A, w) #chart(j,j)
        for l in xrange(2, n+1): 
            for p in xrange(1, n-l+2):
                for t in xrange(1, l):
                    q = p + l - 1
                    d = p + t - 1
                    for A in N:
                        for R in self.binaries.keys():
                            if R[0] == A:
                                B, C = R[1:]
                                if inside[p, B, d] and inside[d+1, C, q]:
                                    inside[p, A, q] = inside[p, A, q] + self.query_x_y_z(A, B, C) * inside[p, B, d] * inside[d+1, C, q]
                    pass
        if inside[1, 'S', n]:
            return inside    
    def recover_tree(self, x, bp, i, j, X):
        if i == j:
            return [X, x[i]]
        else:
            Y, Z, s = bp[i, j, X]
            return [X, self.recover_tree(x, bp, i, s, Y), 
                       self.recover_tree(x, bp, s+1, j, Z)]

if __name__ == "__main__":
    s = 'a boy with a telescope saw a girl'
    parser = PCFGParser('model.bin') 
    parser.writeOutResult(s, '2.out')