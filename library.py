from pickle import NONE


class algcom:
    def __init__(self, n, icod, idet, A, b, tolm):
        self.n, self.icod, self.idet, self.A, self.b, self.tolm = n, idet, icod, A, b, tolm


    def back_substitution(self, y):
        x = [0]*self.n
        x[self.n-1] = y[self.n-1]/self.A[self.n-1][self.n-1]
        for i in range(self.n-1, -1, -1):
            som = 0
            for j in range(i+1, self.n):
                som += x[j]*self.A[i][j]
            x[i] = (y[i] - som)/self.A[i][i]
        return x


    def foward_substitution(self, y):
        x = []
        x.append(y[0]/self.A[0][0])
        for i in range(1, self.n):
            som = 0
            for j in range(0, i):
                som += x[j]*self.A[i][j]
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
    # def solve(self):
    #     if self.icod == 1:
    #         self._LU_decomposition(self)


# n = int(input('Qual o tamanho da matriz?'))
# icod = int(input('Qual o código da operação?'))
# idet = int(input('Devo calcular a determinante?'))
A = eval(input('Insira a matriz A:'))
y = eval(input('Insira o matriz y:'))

# b = eval(input('Insira o vetor b'))
# tolm = int(input('Qual o número máximo de iterações?'))

teste = algcom(3, None, None, A, None, None)

print(teste.back_substitution(y))


