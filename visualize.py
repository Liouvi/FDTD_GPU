import numpy as np
import cv2
import os
import moviepy.video.io.ImageSequenceClip
import matplotlib.pyplot as plt

image_folder=os.getcwd() +"/outputs/"
print(image_folder)
image_files = sorted([os.path.join(image_folder,img)
               for img in os.listdir(image_folder)
               if img.endswith(".dat")])

for i in range(len(image_files)-1):
    Hz = np.loadtxt(f"outputs/Hz{i}.dat")
    Hz = Hz.reshape(1024,1024)
    plt.imsave(f"outputs/"+"Hz" + f"{i}".zfill(5)+ ".png", Hz.T, cmap='RdBu',vmax=1,vmin=-1)

image_files = sorted([os.path.join(image_folder,img)
               for img in os.listdir(image_folder)
               if img.endswith(".png")])

fps=10

clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
clip.write_videofile('outputs/wave.mp4')