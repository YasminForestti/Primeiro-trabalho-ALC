import numpy as np

class algcom:
    def __init__(self, n, icod, idet, A, b, tolm):
        self.n, self.icod, self.idet, self.A, self.b, self.tolm = n, icod, idet, np.array(A), b, tolm
        self.error = 0
        self.autovetor = [0]*self.n
        
    def Jacobi():
        pass
    def PowerMethod(self):
        x = np.array([1]*self.n)
        lamda_ant = 1
        R = np.inf
        maxitter = 1000
        while R > self.tolm or maxitter:
            x  = self.A@x
            lamda = max(x)
            x = x/lamda
            R = abs(lamda - lamda_ant)/abs(lamda)
            lamda_ant = lamda
            maxitter -= 1
        if maxitter == 0: self.error = "Houve divergência no método Power"
        
        return (lamda, x)    
    
    
teste = algcom(3, 3, 0, [[1, 0.2, 0], [0.2, 1, 0.5], [0, 0.5, 1]],[20000, 20000, 20000], 10**(-3))

print(teste.PowerMethod())