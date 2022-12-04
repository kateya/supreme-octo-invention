import matplotlib.pyplot as plt
import matplotlib.animation as animation

def main():
    data = []
    min = 0
    max = 0
    with open("/Users/user/parallel22/results_32_num.txt") as file:
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
    for i in range(32):
        for j in range(32):
            for k in range(32):
                xs.append(i / 32)
                ys.append(j / 32)
                zs.append(k / 32)

    scat = ax.scatter(xs, ys, zs, marker='o', linewidths=0, alpha=0.25, c=data[0], cmap ='viridis', norm=plt.Normalize(min,max))

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    ani = animation.FuncAnimation(fig, update_plot, frames=range(0, numframes), interval=800,
                                  fargs=(data, scat, text))
    #ani.save('./32_num.gif', writer='imagemagick', fps=2)
    plt.show()

def update_plot(i, data, scat, text):
    scat.set_array(data[i])
    text.set_text(f"t = {i}")
    return scat,

main()