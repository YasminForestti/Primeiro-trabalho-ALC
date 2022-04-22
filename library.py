from pickle import NONE
from turtle import forward


class algcom:
    def __init__(self, n, icod, idet, A, b, tolm):
        self.n, self.icod, self.idet, self.A, self.b, self.tolm = n, idet, icod, A, b, tolm


    def _back_substitution(self, y):
        x = [0]*self.n
        x[self.n-1] = y[self.n-1]/self.A[self.n-1][self.n-1]
        for i in range(self.n-1, -1, -1):
            som = 0
            for j in range(i+1, self.n):
                if self.icod == 1:
                    som += x[j]*self.A[i][j]
                else:
                    som += x[j]*self.A[j][i]
            x[i] = (y[i] - som)/self.A[i][i]
        return x

    def Erros():
        pass

    def _foward_substitution(self, y):
        x = []
        if self.icod == 1:
            x.append(y[0])   
        else:
            x.append(y[0]/self.A[0][0])
        for i in range(1, self.n):
            som = 0
            for j in range(0, i):
                som += x[j]*self.A[i][j]
            if self.icod == 1:    
                x.append((y[i] - som))
            else: 
                x.append((y[i] - som)/self.A[i][i])    
        return x
            
    def _LU_decomposition(self):
        for i in range(0, self.n - 1):
            for c in range(i, self.n):
                for r in range(i+1, self.n):
                    if c == i:
                        self.A[r][c] = self.A[r][c]/self.A[c][c]
                    else:
                        self.A[r][c] = self.A[r][c] - self.A[r][i] * self.A[i][c]
        return self.A
    
    def Cholesky(self):
        for i in range(self.n):
            self.A[i][i] = (self.A[i][i] - sum([self.A[i][k]**2  for k in range(0,  i)]))**0.5
            for j in range(i+1, self.n):
                self.A[j][i] = (self.A[i][j] -  sum([self.A[i][k]*self.A[j][k]  for k in range(0,  i)]))/self.A[i][i]
                
        return self.A
    
    def solve(self):
        if self.icod == 1:
            self._LU_decomposition()
            det = 1 
            if self.idet != 0:
                for x in range(0,self.n):
                    det *= self.A[x][x]
            return self._back_substitution(self._foward_substitution(self.b)), det
        if self.icod == 2:
            self.Cholesky()
            det = 1 
            if self.idet != 0:
                for x in range(0,self.n):
                    det *= self.A[x][x]
            return  self._back_substitution(self._foward_substitution(self.b)), det**2                
    def output(self):
        with open("output.txt", "w") as arquivo:
            if self.idet != 0:
               sol = self.solve() 
               arquivo.write(f"Solucao do sistema X: {sol[0]}\n Determinante: {sol[1]}\n")
            else: 
                arquivo.write(f"Solucao do sistema X: {sol[0]}\n")
            

# n = int(input('Qual o tamanho da matriz?'))
# icod = int(input('Qual o código da operação?'))
idet = int(input('Devo calcular a determinante?'))
A = eval(input('Insira a matriz A:'))
#y = eval(input('Insira o matriz y:'))

b = eval(input('Insira o vetor b'))
# tolm = int(input('Qual o número máximo de iterações?'))

teste = algcom(3,2, idet, A,b, None)

teste.output()


