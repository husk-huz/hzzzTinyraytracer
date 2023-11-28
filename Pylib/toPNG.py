from PIL import Image
import numpy as np

def read_image(i):
# 读取像素数据
    with open(f'outfiles/{i}.tga', 'r') as f:
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

for i in range(60):
    read_image(i)
    print(f'{i}.png saved')