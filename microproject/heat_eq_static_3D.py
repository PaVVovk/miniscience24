from fenics import *
import numpy as np
import matplotlib.pyplot as plt

# Импортируем готовую mesh-сетку и строим на ней функциональное пространство
mesh = Mesh('cylinder.xml')
x = SpatialCoordinate(mesh)

V = FunctionSpace(mesh, 'P', 1)
W = FunctionSpace(mesh, "DG", 0) #специальное пространство для задания разрывных величин

dist = sqrt(x[0]*x[0] + x[1]*x[1])

#Параметры задачи
R = 7 * 10**(-3)         #диаметр цилиндра, м
r = 5 * 10**(-4)         #диаметр проволоки, м
L = 0.4                  #длина проволоки, м
T_0 = 23.0               #температура, поддерживаемая термостатом
kappa_wire = 71          #коэффициент теплопроводности проволоки, Вт/(м * К)

#Коэффициент теплопроводности, Вт/(м * К)
kappa = project(conditional(le(dist, r/2), kappa_wire, 0.025), W)  

r_0 = 11 * 10**(-8)      #удельное сопротивление проволоки, Ом*м
alpha = 0.0039           #темп. коэф. сопротивления платины, K^(-1)
U = 0.5                  #разность потенциалов на концах проволоки, В

# Граничное условие в форме Дирихле
T_D = Constant(T_0)

tol = 1E-14
def boundary_D(x, on_boundary):
    return (on_boundary)
#если исключить второе условие, то на боковых границах должно выполняться условие Неймана dT/dn = 0

bc = DirichletBC(V, T_D, boundary_D)

#Индикаторная ф-я, определяющая область проволоки
ind = project(conditional(le(dist, r/2), 1, 0), W) 

#Функции, создающие неоднородность

def f(T):
    '''Функция, определяющая объёмную мощность теплоты, выделяемой в проволоке'''
    return (U / L)**2 / (r_0 * (1 + alpha * (T - T_0)))

def kappa_lin(T):
    '''Линейная часть приращения коэф. теплопроводности воздуха'''
    return  10**(-4) * (T - T_0)

# Постановка вариационной задачи
T = TrialFunction(V) 
v = TestFunction(V)
T_ = interpolate(T_D, V)  #начальное приближение для нелинейной задачи

F = ((kappa + (1 - ind) * kappa_lin(T)) * dot(grad(T), grad(v)) - ind * f(T) * v)*dx
F  = action(F, T_)

J  = derivative(F, T_, T)

problem = NonlinearVariationalProblem(F, T_, bc, J)
solver  = NonlinearVariationalSolver(problem)

res_file = File('heat3D_static/solution.pvd')

# Считаем решение
solver.solve()

# Сохраняем в отдельный файл VTK
res_file << T_

#График распределения температур вдоль проволоки
X = np.linspace(0, L, 1000)
temps = [T_(Point(0, 0, x)) for x in X]
plt.figure()
plt.plot(X * 10**2, temps)
plt.title("График распределения температур вдоль проволоки")
plt.xlabel("x, см")
plt.ylabel("T(x), \u00b0 C")
plt.grid(visible=True)
plt.show()

