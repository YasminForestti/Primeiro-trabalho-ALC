import numpy as np

class algcom:
    def __init__(self, icod, N, pares, x):
        self.icod, self.N, self.pares, self.x = icod, N, pares, x
        self.erro = ""
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
        return sum([phis[i]*self.pares[i][1] for i in range(self.N)])
        


    def Regressao(self):
        P_1 = np.zeros((self.N, 4))
        P_2 = np.zeros((self.N, 2))
        P_3 = np.zeros((self.N, 2))
        Y = np.array([par[1] for par in self.pares])

        f_1 = [lambda x: 1, lambda x: x, lambda x: x**2, lambda x: np.e**x]
        f_2 = [lambda x: 1, lambda x: x]
        f_3 = [lambda x: 1/(np.e**x), lambda x:np.log(x)]
        if self.N >= 4:
            for i in range(self.N):
                for j in range(4):
                    P_1[i][j] = f_1[j](self.pares[i][0])
            for i in range(self.N):
                for j in range(2):
                    P_2[i][j] = f_2[j](self.pares[i][0])
                    if self.pares[i][0] <= 0:
                        self.erro = 3
                    else:
                        P_3[i][j] = f_3[j](self.pares[i][0])
        elif self.N >= 2:
            self.erro = 1
            for i in range(self.N):
                for j in range(2):
                    P_2[i][j] = f_2[j](self.pares[i][0])
                    if self.pares[i][0] <= 0:
                        self.erro = 2
                    else:
                        P_3[i][j] = f_3[j](self.pares[i][0])
        else:
            self.erro = 0
            return []
           
        if self.erro == 1:
            B_2 = np.linalg.inv(P_2.T@P_2)@(P_2.T)@Y
            B_3 = np.linalg.inv(P_3.T@P_3)@(P_3.T)@Y
            y2 = sum([B_2[i] * f_2[i](self.x) for i in range(len(f_2))])
            y3 = sum([B_3[i] * f_3[i](self.x) for i in range(len(f_3))])
            return {"2": y2, "3": y3}
        elif self.erro == 2:
            B_2 = np.linalg.inv(P_2.T@P_2)@(P_2.T)@Y
            y2 = sum([B_2[i] * f_2[i](self.x) for i in range(len(f_2))])
            return {"2": y2}
        elif self.erro == 3:
            B_1 = np.linalg.inv(P_1.T@P_1)@(P_1.T)@Y
            B_2 = np.linalg.inv(P_2.T@P_2)@(P_2.T)@Y
            y1 = sum([B_1[i] * f_1[i](self.x) for i in range(len(f_1))])
            y2 = sum([B_2[i] * f_2[i](self.x) for i in range(len(f_2))])
            return {"1": y1, "2": y2}
        else:
            B_1 = np.linalg.inv(P_1.T@P_1)@(P_1.T)@Y
            B_2 = np.linalg.inv(P_2.T@P_2)@(P_2.T)@Y
            B_3 = np.linalg.inv(P_3.T@P_3)@(P_3.T)@Y
            y1 = sum([B_1[i] * f_1[i](self.x) for i in range(len(f_1))])
            y2 = sum([B_2[i] * f_2[i](self.x) for i in range(len(f_2))])
            y3 = sum([B_3[i] * f_3[i](self.x) for i in range(len(f_3))])
            return  {"1": y1, "2": y2, "3": y3}

    def output(self):
        with open("output.txt", "w") as arquivo:
            if self.icod == 1:
                y = self.Lagrange()
                text = f"f({self.x}) = {y}"
                arquivo.write(text)
            elif self.icod == 2:
                ys = self.Regressao()
                text = ""

                for k, v in ys.items():
                    text += f"f{k}({self.x}) = {v}\n"
                if self.erro == 0:
                    text += "Numero de pontos insuficiente para fazer regressao com o conjunto de funcoes escolhidas"
                elif self.erro == 1:
                    text += "Numero de pontos insuficiente para fazer regressao com o primeiro conjunto de funcoes escolhidas"
                elif self.erro == 2:
                    text +="Numero de pontos insuficiente para fazer regressao com o primeiro conjunto de funcoes escolhidas \n"
                    text +="Nao e possivel fazer regressao com o terceiro conjunto de funcoes escolhidas \n"
                elif self.erro == 3:
                    text +="Nao e possivel fazer regressao com o terceiro conjunto de funcoes escolhidas"
                arquivo.write(text)

    
with open('pares.txt', 'r') as arq:
    text = arq.readlines()
    pares = []
    for line in text:
        par = line.split(' ')
        if par[1][-1] == "\n":
            par[1] = par[1][:-1]
        par = (float(par[0]), float(par[1]))
        pares.append(par)

icod = int(input("Qual o codigo do metodo? "))
N = int(input("Qual o numero de pontos? "))
x = int(input("Qual o valor de x? "))
task3 = algcom(icod, N, pares, x)
task3.output()