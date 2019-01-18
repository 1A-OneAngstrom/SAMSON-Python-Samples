import datetime
import imageio

def create_gif(filenames, duration):
        '''
        function that creates a gif-file from images thanks to imageio package
        '''
        images = []
        for filename in filenames:
                images.append(imageio.imread(filename))
        output_gif_file = 'imageio-%sam.gif' % datetime.datetime.now().strftime('%Y-%M-%d-%H-%M-%S')  # generate name for the gif-file
        output_mp4_file = 'imageio-%sam.mp4' % datetime.datetime.now().strftime('%Y-%M-%d-%H-%M-%S')  # generate name for the mp4-file
        imageio.mimsave(output_gif_file, images, duration=duration)         # save clip to the gif-file
        imageio.mimsave(output_mp4_file, images)                            # save clip to the mp4-file
