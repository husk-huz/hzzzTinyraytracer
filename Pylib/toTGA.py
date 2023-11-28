from PIL import Image

# 打开要转换的PNG图片
image = Image.open("obj/models/rock/rock.png")

# 保存为TGA图片
image.save("obj/models/rock/rock.tga")