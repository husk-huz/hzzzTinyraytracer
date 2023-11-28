# import imageio
# from PIL import Image

# # 读取图片序列，将图片文件名按照数字顺序排序
# array = []

# for i in range(30):
#     array.append("outputs/" + str(i) + ".tga")

# image_files = sorted(array)

# # 读取第一张图片，获取宽度和高度信息
# with Image.open(image_files[0]) as im:
#     width, height = im.size

# # 初始化一个GIF编写器
# with imageio.get_writer('/new_output.gif', mode='I', duration=0.03) as writer :
#     # 写入每一张图片
#     for image_file in image_files:
#         with Image.open(image_file) as im:
#             # 将所有图片的大小调整为第一张图片的大小
#             # im = im.resize((width, height))
#             # 将图片转换为RGB格式
#             im = im.convert('RGB')
#             # 将图片转换为numpy数组
#             im_array = imageio.core.asarray(im)
#             # 写入GIF编写器
#             writer.append_data(im_array)
#             writer.close()

# print("GIF拼接完成！")
# from PIL import Image, ImageFile
# ImageFile.LOAD_TRUNCATED_IMAGES = True

# import cv2

# # 读取TGA文件
# from PIL import Image

# for i in range(60):
# # 打开TGA图像并转换为RGBA模式
#     with open(f'outputs/{i}.tga', 'rb') as f:
#         # img = Image.open(f).convert('RGBA')
#         img = Image.open(f)
#         # img = Image.open(f'outputs/{i}.tga').convert('RGBA')

#         # 获取图片大小
#         width, height = img.size

#         # 创建一个新的RGBA图像对象，背景为透明
#         new_img = Image.new('RGBA', (width, height), (0, 0, 0, 0))

#         # 遍历原始图像像素，并将黑色替换为透明
#         for x in range(width):
#             for y in range(height):
#                 pixel = img.getpixel((x,y))
#                 if pixel[0] == 0 and pixel[1] == 0 and pixel[2] == 0: # 如果是黑色
#                     new_img.putpixel((x,y), (0, 0, 0, 0)) # 设置为透明
#                 else:
#                     new_img.putpixel((x,y), pixel) # 其他颜色不变

# # 保存为PNG格式
#     new_img.save(f'outputs/{i}.png')

# 将图像保存为PNG格式

# for i in range(30):
#     img = cv2.imread(f'outputs/{i}.tga', cv2.IMREAD_UNCHANGED)
#     cv2.imwrite(f'outputs/{i}.jpg', img)

# import imageio.v2 as imageio
# with imageio.get_writer(uri='test.gif', mode='I', fps=30) as writer:
#     for i in range(60):
#         writer.append_data(imageio.imread(f'outputs/{i}.png'))
#         writer.close()
#         imageio.mimsave("hello.gif", gif_images, fps=5)
import imageio.v2 as imageio

gif_images = []
for i in range(120):
    gif_images.append(imageio.imread("outputs/" + str(i) + ".png"))   # 读取多张图片
imageio.mimsave("test.gif", gif_images, fps=45)   # 转化为gif动画


