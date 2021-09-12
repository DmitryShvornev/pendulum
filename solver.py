import json

with open("task.json", 'r') as task_file:
	task = json.load(task_file)
	
print(task)

'''
y_n - решение системы уравнений в текущей точке x_n
порядок следования уравнений такой же, как в бумажке:
theta1, theta2, omega1, omega2. 
функция f_func  возвращает на выход вектор
'''
def f_func(x_n, y_n):
	m1 = task["segments"][0]["m"]
	m2 = task["segments"][1]["m"]
	l1 = task["segments"][0]["l"]
	l2 = task["segments"][1]["l"]

	A  = 9*m2*(cos(y_n[0]-y_n[1]))**2-16*m1-48*m2
	B = 9*m2*cos(y_n[0]-y_n[1])*sin(y_n[0]-y_n[1])
	C = 18*9.81*cos(y_n[0]-y_n[1])
	D = 12*sin(y_n[0]-y_n[1])
	E = 24*9.81*(m1+2*m2)

	omega1_n_1 = (B*l1*(y_n[2])**2 + C*m2*sin(y_n[1]) + D*m2*l2*y_n[3]**2 - E*sin(y_n[0])) / (l1*A)
	omega2_n_1 = (-B*l2*y_n[3]**2 + C*(m1+2*m2)*sin(y_n[0]) - D*(m1+3*m2)*l1*y_n[2]**2 - E*sin(y_n[1]) - 24*9.81*m2*sin(y_n[1]))/(l2*A)

	return [y_n[2],y_n[3],omega1_n_1,omega1_n_2]


def runge_kutta(dt, f1, f2):
	theta1_ = []
	theta2_ = []
	t_ = []
	t_prev = 0
	theta1_prev = task["segments"][0]["theta"]
	theta2_prev = task["segments"][1]["theta"]
	theta1_.append(theta1_prev)
	theta2_.append(theta2_prev)
	t_.append(t_prev)
	'''
	предлагаю переделать рунге-кутта  так, чтобы можно было работать с массивом 4 переменных, а не с каждой переменной в отдельности
	+ по идее, надо как-то сохранять результат предыдущих вычислений, чтобы можно было потом решение построить по точкам
	'''

	for i in range(1, 101):
		t = t_prev + dt
		
		k1_1 = f(t_prev, theta1_prev)
		k2_1 = f(t_prev + dt / 2, theta1_prev * k1_1 / 2)
		k3_1 = f(t_prev + dt / 2, theta1_prev * k2_1 / 2)
		k4_1 = f(t_prev + dt, theta1_prev + dt * k3_1)
		k1_2 = f(t_prev, theta1_prev)
		k2_2 = f(t_prev + dt / 2, theta1_prev * k1_1 / 2)
		k3_2 = f(t_prev + dt / 2, theta1_prev * k2_1 / 2)
		k4_2 = f(t_prev + dt, theta1_prev + dt * k3_1)
		
		theta1 = theta1_prev + dt * (k1_1 + 2 * k2_1 + 2 * k3_1 + k4_1) / 6
		theta2 = theta1_prev + dt * (k1_2 + 2 * k2_2 + 2 * k3_2 + k4_2) / 6
		
		t_prev = t
		theta1_prev = theta1
		theta2_prev = theta2
		
		theta1_.append(theta1_prev)
		theta2_.append(theta2_prev)
		t_.append(t_prev)
		
	return theta1_, theta2_
	
	