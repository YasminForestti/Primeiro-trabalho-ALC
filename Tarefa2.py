import numpy as np

class algcom:
    def __init__(self, n, icod, idet, A, tolm):
        self.n, self.icod, self.idet, self.A, self.tolm = n, icod, idet, np.array(A), tolm
        self.error = 0
        self.autovetor = [0]*self.n
        self.maxiter = 1000
    def Jacobi(self):
        X = np.identity(self.n)
        while(self.maxiter):
            max_value = 0
            i_max, j_max = 0, 0
            for i in range(self.n):
                for j in range(self.n):
                    if max_value < abs(self.A[i][j]) and i != j:
                        max_value = abs(self.A[i][j])
                        i_max, j_max = i, j
            if max_value < self.tolm:
                return X
            phi = 0
            if self.A[i_max][i_max] == self.A[j_max][j_max]:
                phi = np.pi/4
            else:
                phi = np.arctan( (2 * self.A[i_max][j_max])/(self.A[i_max][i_max] - self.A[j_max][j_max]))/2
            P = np.identity(self.n)
            P[i_max][i_max] = np.cos(phi)
            P[j_max][j_max] = np.cos(phi)
            P[i_max][j_max] = -np.sin(phi)
            P[j_max][i_max] = np.sin(phi)
            self.A = (P.T)@self.A@P
            X = X@P
            self.maxiter -= 1
        self.error = "Houve divergência no método Jacobi"
        return X
    def PowerMethod(self):
        x = np.array([1]*self.n)
        lamda_ant = 1
        R = np.inf
        while R > self.tolm and self.maxiter:
            x  = self.A@x
            lamda = max(x)
            x = x/lamda
            R = abs(lamda - lamda_ant)/abs(lamda)
            lamda_ant = lamda
            self.maxiter -= 1
        if self.maxiter == 0: 
            self.error = "Houve divergência no método Power"
        return (lamda, x)  
      
    def output(self):
        with open("output.txt", "w") as arquivo:
            if self.icod == 1:
                lamda, x = self.PowerMethod()
                if self.error:
                    arquivo.write(f"Error: {self.error}\n")
                else:
                    text = f"lambda = {lamda} \n X = {x}"
                    arquivo.write(text)
            elif self.icod == 2:
                X = self.Jacobi()
                if self.error:
                    arquivo.write(f"Error: {self.error}\n")
                else:
                    if self.idet == 0:
                        text = f"autovalores na diagonal de: {self.A} \n autovetores nas colunas de: {X}"
                        arquivo.write(text)
                    else:
                        det = 1
                        for i in range(self.n):
                            det *= self.A[i][i]
                        text = f"autovalores na diagonal de: {self.A} \n autovetores nas colunas de: {X} \n O determinante de A: {det}"
                        arquivo.write(text)

                
n = int(input('Qual a ordem da matriz A?'))
icod = int(input('Qual o código da operação?'))
idet = int(input('Devo calcular a determinante?'))
tolm = float(input('Qual o número máximo de iterações?'))

with open('A.txt', 'r') as arq:
    contlinhas = 0
    text = arq.readlines()
    A =  []
    for line in text:
        linha = line.split(' ')
        if linha[n-1][-1] == '\n':
            linha[n-1] = linha[n-1][:-1]
        for i in range(len(linha)):
            linha[i] = int(linha[i])
        A.insert(contlinhas,linha)    
        contlinhas+=1
       
teste = algcom(n, icod, idet, A,tolm)

teste.output()