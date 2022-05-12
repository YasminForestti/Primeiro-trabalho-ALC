
from joblib import PrintTime
import numpy as np

class algcom:
    def __init__(self, icod, N, pares, x):
        self.icod, self.N, self.pares, self.x = icod, N, pares, x
    def Lagrange(self):
        phis = []
        for i in range(self.N):
            prod_num = 1
            prod_den = 1
            for k in range(self.N):
                if k != i:
                    prod_num *= self.x - self.pares[k][0]
                    prod_den *= self.pares[i][0] - self.pares[k][0]
            phis.append(prod_num/prod_den)
        print(phis)
        return sum([phis[i]*self.pares[i][1] for i in range(self.N)])
        


    def Regressao(self):
        P_1 = np.zeros((self.N, 4))
        P_2 = np.zeros((self.N, 2))
        P_3 = np.zeros((self.N, 2))
        Y = np.array([par[1] for par in self.pares])
        f_1 = [lambda x: 1, lambda x: x, lambda x: x**2, lambda x: np.e**x]
        f_2 = [1, lambda x: x]
        f_3 = [lambda x: 1/(np.e**x), lambda x:np.log(x)]
        for i in range(self.N):
            for j in range(4):
                P_1[i][j] = f_1[j](self.pares[0])
        
        for i in range(self.N):
            for j in range(2):
                P_2[i][j] = f_2[j](self.pares[0])
                P_3[i][j] = f_3[j](self.pares[0])
        B_1 = np.linalg.inv(P_1.T@P_1)@(P_1.T)@Y
        B_2 = np.linalg.inv(P_2.T@P_2)@(P_2.T)@Y
        B_3 = np.linalg.inv(P_3.T@P_3)@(P_3.T)@Y

        return B_1, B_2, B_3

    def output(self):
        with open("output.txt", "w") as arquivo:
            if self.icod == 1:
                y = self.Lagrange()
                text = f"f({self.x}) = {y}"
                arquivo.write(text)
            elif self.icod == 2:
                pass
                

    
with open('tarefa3.txt', 'r') as arq:
    text = arq.readlines()
    pares = []
    for line in text:
        par = line.split(', ')
        if par[1][-1] == "\n":
            par[1] = par[1][:-1]
        par = (int(par[0]), int(par[1]))
        pares.append(par)



teste = algcom(1, 3, pares, 2)

teste.output()