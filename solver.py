import json

with open("task.json", 'r') as task_file:
	task = json.load(task_file)
	
print(task)

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
	
	