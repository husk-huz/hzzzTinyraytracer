from PIL import Image
import numpy as np
import imageio.v2 as imageio


numbers = 60

def read_image(i):
# 读取像素数据
    with open(f'outfiles/{i}.hzzz', 'r') as f:
        pixels = f.read().splitlines()
        # print(pixels[0])

    # 转换数据格式
    pixels = [tuple(map(int, p.split())) for p in pixels]
    width = pixels[0][0]
    height = pixels[0][1]

    pixels.remove(pixels[0])

    # pixels = np.array(pixels).reshape(height, width, 3)
    npimage = np.zeros((height, width, 3), dtype=np.uint8)

    for x in range(width):
        for y in range(height):
            npimage[x][y][0] = pixels[y][3*x  ]
            npimage[x][y][1] = pixels[y][3*x+1]
            npimage[x][y][2] = pixels[y][3*x+2]

    # 创建Image对象
    # img = Image.new('RGB', (width, height), 'white')
    img = Image.fromarray(np.flip(npimage, 0), 'RGB')
    # img.putdata(pixels)

    # 保存为png图片
    img.save(f'outputs/{i}.png', 'PNG')

def read_and_save_image(numbers):
    gif_images = []
    for i in range(numbers):
        with open(f'outfiles/{i}.hzzz', 'r') as f:
            pixels = f.read().splitlines()
            # print(pixels[0])

        # 转换数据格式
        pixels = [tuple(map(int, p.split())) for p in pixels]
        width = pixels[0][0]
        height = pixels[0][1]

        pixels.remove(pixels[0])

        # pixels = np.array(pixels).reshape(height, width, 3)
        npimage = np.zeros((height, width, 3), dtype=np.uint8)

        for x in range(width):
            for y in range(height):
                npimage[x][y][0] = pixels[y][3*x  ]
                npimage[x][y][1] = pixels[y][3*x+1]
                npimage[x][y][2] = pixels[y][3*x+2]

        # 创建Image对象
        # img = Image.new('RGB', (width, height), 'white')
        img = Image.fromarray(np.flip(npimage, 0), 'RGB')
        gif_images.append(img)   
        print("image " + str(i) + " read")
    imageio.mimsave("test.gif", gif_images, duration=20)   
    print("test.gif saved")

read_and_save_image(numbers)

# for i in range(numbers):
#     read_image(i)
#     print(f'{i}.png saved')

# gif_images = []
# for i in range(numbers):
#     gif_images.append(imageio.imread("outputs/" + str(i) + ".png"))   # 读取多张图片
# imageio.mimsave("test.gif", gif_images, duration=20)   # 转化为gif动画

# print("test.gif saved")

