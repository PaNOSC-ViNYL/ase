from ase.visualize.plot import animate

def write_gif(filename, images, interval=200, save_count=100,
              show_unit_cell=2, **parameters):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    animation = animate(images, fig=fig, ax=ax,
                        interval=interval, save_count=save_count,
                        show_unit_cell=show_unit_cell,
                        **parameters)
    #animation.save(filename, writer='imagemagick')
    #animation.save(filename + '.mp4')
