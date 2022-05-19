class algcom:
    def __init__(self, n, icod, idet, A, b, tolm):
        self.n, self.icod, self.idet, self.A, self.b, self.tolm = n, icod, idet, A, b, tolm
        self.error = 0

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
                        #primeira iteração
                        try:
                            self.A[r][c] = self.A[r][c]/self.A[c][c] # definindo alpha e betha positivos
                        except:
                            # Acontece quando o pivô é zero
                            self.error = "Decomposição LU não é possível"
                            return
                    else:
                        self.A[r][c] = self.A[r][c] - self.A[r][i] * self.A[i][c]
                        
        return
    
    def Cholesky(self):
        for i in range(self.n):
            if self.A[i][i] <= 0:
                self.error = "Esta matriz nao e positiva definida, nao e possivel aplicar Cholesky"
                return
            for j in range(self.n):
                if self.A[i][j] != self.A[j][i]:
                    self.error = "Esta matriz nao e simetrica, nao e possivel aplicar Cholesky"
                    return
        for i in range(self.n):
            raiz = self.A[i][i] - sum([self.A[i][k]**2  for k in range(0,  i)])
            if raiz < 0:
                # podemos n ter detectado que a matriz não é positiva definida
                self.error =  "Decomposição Cholesky não é possível"
                return
            else:
                self.A[i][i] = (raiz)**0.5
            for j in range(i+1, self.n):
                #preenchendo a inferior e consultando a superior
                self.A[j][i] = (self.A[i][j] -  sum([self.A[i][k]*self.A[j][k]  for k in range(0,  i)]))/self.A[i][i]

        return
    

    def norma(self, x, y):
        return (sum([(x[i] - y[i])**2 for i in range(len(x))])**0.5)/(sum([x[i]**2 for i in range(len(x))])**0.5)


    def Jacobi(self):
        x = [1]*self.n
        maxitter = 1000
        temp = x.copy()
        while(maxitter):
            maxitter -= 1
            for i in range(self.n):
                soma = 0
                for j in range(self.n):
                    if j != i:
                        soma += self.A[i][j] * x[j]
                temp[i] = (self.b[i] - soma)/self.A[i][i]

            res = self.norma(temp, x)
            x = temp.copy()
            if res <= self.tolm:
                return x
        self.error = "Houve divergência no método Jacobi"

    def GaussSaiden(self):
        x = [1]*self.n
        maxitter = 1000
        
        while(maxitter):
            x_ant = x.copy()
            maxitter -= 1
            for i in range(self.n):
                soma = 0
                for j in range(self.n):
                    if j != i:
                        soma += self.A[i][j] * x[j]
                x[i] = (self.b[i] - soma)/self.A[i][i]

            res = self.norma(x,x_ant)
            if res <= self.tolm:
                
                return x
        self.error = "Houve divergência no método Jacobi"

    def solve(self):
        if self.icod == 1:
            self._LU_decomposition()
            det = 1 
            if self.idet != 0:
                for x in range(0,self.n):
                    det *= self.A[x][x]
            if not self.error:
                # LY = B
                # UX = Y
                return self._back_substitution(self._foward_substitution(self.b)), det

        if self.icod == 2:
            self.Cholesky()
            det = 1 
            if self.idet != 0:
                for x in range(0,self.n):
                    det *= self.A[x][x]
            if not self.error:
                # LY = B
                # UX = Y
                # U = L.T
                return self._back_substitution(self._foward_substitution(self.b)), det**2
        if self.icod == 3:
            if not self.error:
                return self.Jacobi(), 1
    def output(self):
        sol = self.solve() 
        with open("output.txt", "w") as arquivo:
            if not self.error:
                if self.idet != 0:
                    arquivo.write(f"Solucao do sistema X: {sol[0]}\n Determinante: {sol[1]}\n")
                else: 
                    arquivo.write(f"Solucao do sistema X: {sol[0]}\n")
            else:
                 arquivo.write(f"Error: {self.error}\n")

n = int(input('Qual o tamanho da matriz?'))
icod = int(input('Qual o código da operação?'))
idet = int(input('Devo calcular a determinante?'))
tolm = float(input('Qual a tolerância?'))

with open('A.txt', 'r') as arq:
    contlinhas = 0
    text = arq.readlines()
    A =  []
    for line in text:
        linha = line.split(' ')
        if linha[n-1][-1] == '\n':
            linha[n-1] = linha[n-1][:-1]
        for i in range(len(linha)):
            linha[i] = float(linha[i])
        A.insert(contlinhas,linha)    
        contlinhas+=1
    
    with open('b.txt', 'r') as arq:
        contlinhas = 0
        text = arq.readlines()
        b =  []
        for line in text:
            linha = line.split(' ')
            if linha[0][-1] == '\n':
                linha[0] = linha[0][:-1]            
            b.insert(contlinhas,float(linha[0]))    
            contlinhas+=1
        
task1 = algcom(n, icod, idet, A, b,tolm)

task1.output()