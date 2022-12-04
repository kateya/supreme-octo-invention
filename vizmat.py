import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

n = 2

def main():
    data = []
    min = 0
    max = 0
    with open("./results_1_32.txt") as file:
        for line in file:
            iter = []
            for num in line.split():
                num = float(num)
                if num > max:
                    max = num
                if num < min:
                    min = num
                iter.append(num)
            data.append(iter)

    numframes = len(data)

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(projection='3d')
    text = ax.set_title('t = 0')

    xs = []
    ys = []
    zs = []
    for i in range(33):
        for j in range(33):
            for k in range(33):
                xs.append(i)
                ys.append(j)
                zs.append(k)

    scat = ax.scatter(xs, ys, zs, marker='o', linewidths=0, alpha=0.25, c=data[0], cmap ='viridis', norm=plt.Normalize(min,max))

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    ani = animation.FuncAnimation(fig, update_plot, frames=range(1, numframes), interval=1000,
                                  fargs=(data, scat, text))
    #ani.save('./animation.gif', writer='imagemagick', fps=30)
    plt.show()

def update_plot(i, data, scat, text):
    scat.set_array(data[i])
    text.set_text(f"t = {i}")
    return scat,

main()