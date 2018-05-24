## this goes on the todo list

import pandas as pd
from maxijdefs import *
from matplotlib import pyplot as plt

def checkplot(night, fIDs, path='/media/amanda/demeter/maxi_j1820_070/'):
plt.register_cmap(name='SkyCMAP', data=cdict4)

# if int(id) % savecheck == 0:
#     ## plot tycho star with an aperture drawn on for each image
#     fig = plt.figure(facecolor = 'w')
#     ax = plt.gca()
#     ax.cla()  # clear things for fresh plot
#     plt.imshow(np.arcsinh(im0+10000),cmap='gray')


for i in refs.index:
    # for ref star positions rx, ry (with tycho star gaussian shifts):
    rx = np.rint(refs.loc[i]['X(FITS)']) + dx
    ry = np.rint(refs.loc[i]['Y(FITS)']) + dy

    sky[i], flux[i], skyplot = measure_star(im0, rx, ry, sg)

    # plot sky apertures:
    # plt.imshow(skyplot, cmap='SkyCMAP')

    # to plot apertures:
    # if int(id) % savecheck == 0:
    #     if i == 1:
    #         # maxij in blue
    #         ax.add_artist(plt.Circle((rx, ry), 3. * sg, color='b', fill=False))
    #         plt.text(rx, ry+10*sg, 'maxi j1820+070', ha="center", family='sans-serif', size=14, color = 'blue')
    #         plt.text(10, 235, 'BH f=' + str(np.rint(flux[i])) + ', sky=' + str(np.rint(sky[i])), family='sans-serif', size=14, color='blue')
    #     else:
    #         # the reference stars in red
    #         ax.add_artist(plt.Circle((rx, ry), 3. * sg, color='r', fill=False))
    #         plt.text(rx, ry-5*sg, 'r' + str(i+1), ha="center", family='sans-serif', size=14, color = 'red')

    ## save checkimages with filenames:
    # if int(id)%savecheck == 0:
    #     # plt.show()
    #     plt.text(10, 220, dataf.loc[id]['filename'], family='sans-serif', size=14, color='yellow')
    #     fig.savefig(pathnam + 'checkims/checkapertures_' + id + '.png')
