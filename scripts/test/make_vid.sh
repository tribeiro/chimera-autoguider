ffmpeg -framerate 1 -i img_%04d.png -c:v libx264 -vf fps=25 -pix_fmt yuv420p guider.mp4
