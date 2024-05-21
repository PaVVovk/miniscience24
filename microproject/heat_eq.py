from fenics import *
import numpy as np
import matplotlib.pyplot as plt

time = 20                  #время процесса
num_steps = 30 * time      #число шагов
dt = time / num_steps      #размер одного шага

# Импортируем готовую mesh-сетку и строим на ней функциональное пространство
mesh = Mesh('rect.xml')
x = SpatialCoordinate(mesh)

V = FunctionSpace(mesh, 'P', 1)
W = FunctionSpace(mesh, "DG", 0) #специальное пространство для задания разрывных величин

dist = abs(x[1])

#Параметры задачи
R = 7 * 10**(-3)         #диаметр цилиндра, м
r = 5 * 10**(-4)         #диаметр проволоки, м
L = 0.4                  #длина проволоки, м
T_0 = 23.0               #температура, поддерживаемая термостатом

c = 140                  #удельная теплоёмкость проволоки, Дж/(кг * К)
rho = 21450              #плотность платины, кг/м^3
kappa_wire = 71          #коэффициент теплопроводности проволоки, Вт/(м * К)

#Коэффициент при dT/dt, Дж/(м^3 * К)
coef_dt = project(conditional(le(dist, r/2), c*rho, 2.5 * 10**5 / (273.15 + T_0)), W) 
#Коэффициент теплопроводности, Вт/(м * К)
kappa = project(conditional(le(dist, r/2), kappa_wire, 0.025), W) 

r_0 = 11 * 10**(-8)      #удельное сопротивление проволоки, Ом*м
alpha = 0.0039           #темп. коэф. сопротивления платины, K^(-1)
U = 0.5                  #разность потенциалов на концах проволоки, В

# Граничное условие в форме Дирихле
T_D = Constant(T_0)

tol = 1E-14
def boundary_D(x, on_boundary):
    return (on_boundary 
            and ((near(x[1], -R/2, tol) or near(x[1], R/2, tol))
                 or (near(x[0], 0, tol) or near(x[0], L, tol))))
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

# Начальные значения величины T
T_n = interpolate(T_D, V) #интерполируем граничное условие на всё V

# Постановка вариационной задачи
T = TrialFunction(V) 
v = TestFunction(V)
T_ = interpolate(T_D, V)  #начальное приближение для нелинейной задачи

F = ((kappa + (1 - ind) * kappa_lin(T))*dt*dot(grad(T), grad(v)) + (coef_dt*(T - T_n) - ind * f(T)*dt)*v)*x[1]*dx
F  = action(F, T_)

J  = derivative(F, T_, T)

problem = NonlinearVariationalProblem(F, T_, bc, J)
solver  = NonlinearVariationalSolver(problem)

res_file = File('heat2D/solution.pvd')

# Делаем шаги по времени
#t = 0; здесь переменная времени явно нигде не используется
for n in range(num_steps):

    # Считаем решение для отдельного шага
    solver.solve()

    # Сохраняем в отдельный файл VTK
    res_file << T_

    # Обновляем предыдущее решение
    T_n.assign(T_)

#График распределения температур вдоль проволоки
X = np.linspace(0, L, 1000)
temps = [T_(Point(x, 0)) for x in X]
plt.figure()
plt.plot(X * 10**2, temps)
plt.title(f"График распределения температур вдоль проволоки, t = {time} с")
plt.xlabel("x, см")
plt.ylabel("T(x), \u00b0 C")
plt.grid(visible=True)
plt.show()
