import json
from math import *
import matplotlib.pyplot as plt

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

	A  = 9 * m2 * (cos(y_n[0] - y_n[1])) ** 2 - 16 * m1 - 48 * m2
	B = 9 * m2 * cos(y_n[0] - y_n[1]) * sin(y_n[0] - y_n[1])
	C = 18 * 9.81 * cos(y_n[0] - y_n[1])
	D = 12 * sin(y_n[0] - y_n[1])
	E = 24 * 9.81 * (m1 + 2 * m2)

	omega1_n_1 = (B * l1 * (y_n[2]) ** 2 + C * m2 * sin(y_n[1]) + D * m2 * l2 * y_n[3] ** 2 - E * sin(y_n[0])) / (l1 * A)
	omega2_n_1 = (-B * l2 * y_n[3] ** 2 + C * (m1 + 2 * m2) * sin(y_n[0]) - D * (m1 + 3 * m2) * l1 * y_n[2] ** 2 - E * sin(y_n[1]) - 24 * 9.81 * m2 * sin(y_n[1])) / (l2 * A)

	return [y_n[2],y_n[3],omega1_n_1,omega2_n_1]


def runge_kutta(dt, f):
	theta1_ = []
	theta2_ = []
	omega1_ = []
	omega2_ = []
	t_ = []
	t_prev = 0
	theta1_prev = task["segments"][0]["theta"]
	theta2_prev = task["segments"][1]["theta"]
	omega1_prev = task["segments"][0]["omega"]
	omega2_prev = task["segments"][1]["omega"]
	theta1_.append(theta1_prev)
	theta2_.append(theta2_prev)
	omega1_.append(omega1_prev)
	omega2_.append(omega2_prev)
	t_.append(t_prev)

	for i in range(100):
		t = t_prev + dt
		
		
		k1 = f(t_prev, [theta1_[i], theta2_[i], omega1_[i], omega2_[i]])
		k2 = f(t_prev + dt / 2, [theta1_[i] * k1[0] / 2, theta2_[i] * k1[1] / 2, omega1_[i] * k1[2] / 2, omega2_[i] * k1[3] / 2])
		k3 = f(t_prev + dt / 2, [theta1_[i] * k2[0] / 2, theta2_[i] * k2[1] / 2, omega1_[i] * k2[2] / 2, omega2_[i] * k2[3] / 2])
		k4 = f(t_prev + dt, [theta1_[i] + dt * k3[0], theta2_[i] + dt * k3[1], omega1_[i] + dt * k3[2], omega2_[i] + dt * k3[3]])
		
		theta1 = theta1_prev + dt * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6
		theta2 = theta2_prev + dt * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6
		omega1 = omega1_prev + dt * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6
		omega2 = omega2_prev + dt * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]) / 6
		
		t_prev = t
		theta1_prev = theta1
		theta2_prev = theta2
		omega1_prev = omega1
		omega2_prev = omega2
		
		theta1_.append(theta1_prev)
		theta2_.append(theta2_prev)
		omega1_.append(omega1_prev)
		omega2_.append(omega2_prev)
		t_.append(t_prev)
		
	return theta1_, theta2_, omega1_, omega2_, t_


output_data = runge_kutta(0.01, f_func)

fig, ax = plt.subplots()

ax.scatter(output_data[4], output_data[0], c = 'red', label = "theta1")
ax.scatter(output_data[4], output_data[1], c = 'green', label = "theta2")
ax.scatter(output_data[4], output_data[2], c = 'blue', label = "omega1")
ax.scatter(output_data[4], output_data[3], c = 'black', label = "omega2")


ax.set_title('Результаты расчетов')

ax.legend()

fig.set_figwidth(8)
fig.set_figheight(8)  

plt.show()
	
	