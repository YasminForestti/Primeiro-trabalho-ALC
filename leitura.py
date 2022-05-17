from turtle import pensize


n = 3
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

print(A)