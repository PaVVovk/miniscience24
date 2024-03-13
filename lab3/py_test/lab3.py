from fenics import *

T = 5.0            # время процесса
num_steps = 25     # число шагов
dt = T / num_steps # размер одного шага

# Импортируем готовую mesh-сетку и строим на ней функциональное пространство
mesh = Mesh('cat.xml')
V = FunctionSpace(mesh, 'P', 1)

# Граничное условие в форме Дирихле
u_normal = 36.6
u_D = Constant(u_normal)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Начальные значения величины u
u_n = interpolate(u_D, V) #интерполируем граничное условие на всё V

# Постановка вариационной задачи
u = TrialFunction(V)
v = TestFunction(V)

# Создадим источник тепла-"вспышку"
u_hot = u_normal * 5000
x_n = 0.0
y_n = 11.0
z_n = 20.0
sigma = 4.0
f = Expression('A * exp(-1 * (pow(x[0]-x_n, 2) + pow(x[1]-y_n, 2) + pow(x[2]-z_n, 2)) / (2*s*s)) / pow(2*pi*s*s, 1.5)', 
               degree=3, A=u_hot, s=sigma, x_n=x_n, y_n=y_n, z_n=z_n)

F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
a, L = lhs(F), rhs(F)

res_file = File('heat/solution.pvd')

# Делаем шаги по времени
u = Function(V)
t = 0
for n in range(num_steps):

    # Обновляем время (и при необходимости - параметры функций)
    t += dt

    # Считаем решение для отдельного шага
    solve(a == L, u, bc)

    # Сохраняем в отдельный файл VTK
    res_file << u

    # Обновляем предыдущее решение
    u_n.assign(u)