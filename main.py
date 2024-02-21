import tkinter as tk
from tkinter import messagebox
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
import numpy as np
from scipy.integrate import quad
h_eV = 4.136e-15
Const_Planck_Js = 6.626e-34
const_velo = 2.187e6
Velo_Luz = 3e8
Massa_eletron = 9.11e-31
constante_ry = 1.097e7
Massa_proton = 1.673e-27


def menu():
  print('MENU')
  print('''
  1 - Dados sobre FUNÇÃO DE ONDA
  2 - Informações entre niveis
  3 - Probabilidade
  4 - Questões quem envolvem a função de onda 
  5 - Simulação de transição de estado''')

#Determine a função de onda do nível n = X do elétron.
def funçãoOnda():
  nivel = int(input('Digite o nível: '))
  largura = float(
      input('Digite a lagura do poço de potencial infinito em [NM]:'))
  largura1 = largura / 1e9
  psi = math.sqrt(2 / largura1)
  Knivel = (nivel * math.pi) / (largura1)

  print(f'O valor para PSI é: {psi:.4e}')
  print(f'O valor no nivel {nivel} é: {Knivel:.4e}')


def energiaEntreNiveis():
  op = int(input('''[1] Elétron 
[2] Próton
Escolha uma opção: '''))

  if op == 1:

    nivel1 = float(input('Digite o valor para o nivel: '))
    nivel2 = float(input('Digite o valor para o nivel: '))
    largura = float(
        input("Digite a lagura do poço de potencial infinito em [NM]: "))
    largura1 = largura * 1e-9

    Enivel1 = (nivel1**2) * (Const_Planck_Js**
                             2) / (8 * Massa_eletron * largura1**2)
    Enivel2 = (nivel2**2) * (Const_Planck_Js**
                             2) / (8 * Massa_eletron * largura1**2)
    e1j = Enivel1 / 1.602e-19
    e2j = Enivel2 / 1.602e-19
    print('\n')
    print(
        f'Energia para o eletron no nivel {nivel1} = {Enivel1:.4e} tem seu J = {e1j:.4e} eV'
    )
    print("\n")
    print(
        f'Energia para o eletron no nivel {nivel2} = {Enivel2:.4e} tem seu J = {e2j:.4e} eV'
    )
    aux = abs(e1j - e2j)
    print("\n")
    print(f"A energia do fóton absorvido é [eV]: {aux:.4e}")
    print("\n")
    comprimento = (h_eV * Velo_Luz) / aux
    frequencia = (Velo_Luz / comprimento)
    print(f'O comprimento do fóton é [M]: {comprimento:.4e}')
    print("\n")
    print(f"A frequencia do fóton absorvido é [Hz]: {frequencia:.4e}")
    print("\n")
    
    vi = math.sqrt(2 * Enivel1 / Massa_eletron)
    vf = math.sqrt(2 * Enivel2 / Massa_eletron)

    print(f'A velocidade de {Enivel1} é: {vi:.4e} m/s')
    print('\n')
    print(f'A velocidade de {Enivel2} é: {vf:.4e} m/s')
    print("\n")

    comp_i = Const_Planck_Js / math.sqrt(2 * Massa_eletron * Enivel1)
    comp_f = Const_Planck_Js / math.sqrt(2 * Massa_eletron * Enivel2)

    print(f'O comprimento de onda de De Broglie na párticula em {nivel1} é: {comp_i:.4e} m')
    print('\n')
    print(f'O comprimento de onda de De Broglie na párticula em {nivel2} é: {comp_f:.4e} m')
    print('\n')
  
  elif op == 2:

    nivel1 = float(input('Digite o valor para o nivel: '))
    nivel2 = float(input('Digite o valor para o nivel: '))
    print('\n')
    largura = float(
        input("Digite a lagura do poço de potencial infinito em [NM]: "))
    largura1 = largura * 1e-9

    Enivel1 = (nivel1**2) * (Const_Planck_Js**2) / (8 * Massa_proton *
                                                    (largura1**2))
    Enivel2 = (nivel2**2) * (Const_Planck_Js**2) / (8 * Massa_proton *
                                                    (largura1**2))
    e1jp = Enivel1 / 1.602e-19
    e2jp = Enivel2 / 1.602e-19
    print(
        f'Energia para o proton no nivel {nivel1} = {Enivel1:.4e} tem seu J = {e1jp:.4e} eV'
    )
    print("\n")
    print(
        f'Energia para o eletron no nivel {nivel2} = {Enivel2:.4e} tem seu J = {e2jp:.4e} eV'
    )
    aux = abs(e1jp - e2jp)
    print("\n")
    print(f"A energia do fóton absorvido é [eV]: {aux:.4e}")
    print("\n")
    comprimento = (h_eV * Velo_Luz) / aux
    frequencia = (Velo_Luz / comprimento)
    print(f'O comprimento do fóton é [M]: {comprimento:.4e}')
    print("\n")
    print(f"A frequencia do fóton absorvido é [Hz]: {frequencia:.4e}")
    print("\n")

    vi = math.sqrt(2 * Enivel1 / Massa_proton)
    vf = math.sqrt(2 * Enivel2 / Massa_proton)

    print(f'A velocidade de {Enivel1} é: {vi:.3g} m/s')
    print("\n")
    print(f'A velocidade de {Enivel2} é: {vf:.3g} m/s')
    print("\n")

    comp_i = Const_Planck_Js / math.sqrt(2 * Massa_proton * Enivel1)
    comp_f = Const_Planck_Js / math.sqrt(2 * Massa_proton * Enivel2)

    print(f'O comprimento de onda de De Broglie na párticula em {nivel1} é: {comp_i:.4e} m')
    print('\n')
    print(f'O comprimento de onda de De Broglie na párticula em {nivel2} é: {comp_f:.4e} m')
    print('\n')
  else:
    print("Entrada inválida")




def onda(x,L,n):  
  return np.sqrt(2/L) * np.sin(n* np.pi* x/L)
def densidade(x,L,n):
  j = onda(x,L,n)
  return j**2
def probabilidade():
  n = int(input("Qual o nivel que a particula se encontra? "))


  Largura = float(input("Qual a largura do poço? "))
  p = int(input("Qual valor da potência: "))
  print("\n")
  potencia = math.pow(10,p)
  L = Largura*potencia


  minimo = float(input("Qual o nivel minimo? "))
  no= int(input("Qual valor da potência: "))
  print("\n")
  potenci = math.pow(10,no)
  n_min = minimo * potenci


  maximo = float(input("Qual o nivel maximo? "))
  nn = int(input("Qual valor da potência: "))
  print("\n")
  potenci = math.pow(10,nn)
  n_max = maximo * potenci

  probabilidade, _ = quad(densidade, n_min, n_max, args=(L, n))

  porcentagem = probabilidade * 100

  print("A probabilidade de encontrar a partícula é de:", porcentagem)

def larguraENivel():
  psi =  float(input("Entre com o valor de PSI: "))
  largura = (2/(psi**2))*1e9
  print(f'A largura de onda é: {largura:.4e}')
  sin = float(input("Entre com o valor que acompanha o SEN: "))
  nivel = (sin*largura*1e-9)/math.pi
  n = round(nivel, 1)
  print(f"O nivel em que a particula se encontra é: {n}")
  entrada = float(input("Entre com a posição: "))
  conta = (psi**2)*math.sin(sin*entrada*largura*1e-9)**2
  print(f"A probabilidade de encontra-lo na posição {entrada} é: {conta:.4e}")
  n = float(input("Entre com o ESTADO EXCITADO que a particula se encontra: "))
  velo = Const_Planck_Js/(Massa_eletron*largura*1e-9)
  print(f"A velocidade é [M/s]: {velo:.4e}")
  n1 = float(input('Digite o mais alto nivel da transição: '))
  n2 = float(input('Digite o mais baixo nivel da transição: '))
  contaN1= (n1**2)*(Const_Planck_Js**2)/(8*((largura*1e-9)**2)*Massa_eletron)
  contaN2= (n2**2)*(Const_Planck_Js**2)/(8*((largura*1e-9)**2)*Massa_eletron)
  contaN = contaN1 - contaN2
  compri = (Const_Planck_Js*Velo_Luz)/(contaN*1e9)
  print(f'O comprimento é [Nm] {compri:.4e}')

def escolha():
  print('\n')
  A = int(input('Escolha uma opção: '))
  print('\n')
  return A



def calcular():
  try:
      psi = float(psi_entry.get())
      largura = (2 / (psi ** 2)) * 1e9
      largura_label.config(text=f'A largura é: {largura:.4e}')

      sin = float(sin_entry.get())
      nivel = (sin * largura * 1e-9) / math.pi
      n = round(nivel, 1)
      nivel_label.config(text=f"O nível em que a partícula se encontra é: {n}")

      animar_particula(int(n), largura)  

  except ValueError:
      messagebox.showerror("Erro", "Por favor, insira valores válidos.")

def animar_particula(nivel, largura):
  levels = list(range(1, 9))
  photons_emitted = [0] * 8 

  fig, ax = plt.subplots()
  ax.set_xlabel('Largura')
  ax.set_ylabel('Níveis de energia')
  ax.set_title('Animação da Partícula')

  ax.set_xlim(0, largura)  
  ax.set_ylim(0, 9)  

  def animate(frame):
      ax.clear()
      ax.scatter([largura] * 8, levels, s=100, color='lightgray')  # Plotando níveis de energia
      ax.scatter(largura, levels[frame], s=100, color='blue')  # Partícula no nível atual
      ax.set_xlabel('Largura')
      ax.set_ylabel('Níveis de energia')
      ax.set_title(f'Animação da Partícula - Nível {frame + 1}/{nivel}')
      ax.set_xlim(0, largura)  # Mantém o limite do eixo x atualizado
      ax.set_ylim(0, 9)  # Mantém o limite do eixo y atualizado

  ani = animation.FuncAnimation(fig, animate, frames=nivel, blit=False, interval=500)
  plt.show()




while True:

  menu()

  x = escolha()

  if (x == 1):
    funçãoOnda()
    

  if (x == 2):
    energiaEntreNiveis()

  if (x==3):
    probabilidade()

  if (x==4):
      larguraENivel()

  if (x==5):
    root = tk.Tk()
    root.title("Simulação de Saltos Quânticos")

    # Criando os widgets
    psi_label = tk.Label(root, text="Primeiro Valor:")
    psi_entry = tk.Entry(root)
    psi_label.pack()
    psi_entry.pack()

    sin_label = tk.Label(root, text="Valor que acompanha o SEN:")
    sin_entry = tk.Entry(root)
    sin_label.pack()
    sin_entry.pack()

    calcular_button = tk.Button(root, text="Calcular", command=calcular)
    calcular_button.pack()

    largura_label = tk.Label(root, text="")
    largura_label.pack()

    nivel_label = tk.Label(root, text="")
    nivel_label.pack()

    root.mainloop()

  if (x == 0):
    break
