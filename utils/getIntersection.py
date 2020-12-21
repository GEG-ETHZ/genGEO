# Licensed under LGPL 2.1, please see LICENSE for details
# https://www.gnu.org/licenses/lgpl-2.1.html
#
# The work on this project has been performed at the GEG Group at ETH Zurich:
# --> https://geg.ethz.ch
#
# The initial version of this file has been implemented by:
#
#     Philipp Schaedle (https://github.com/philippschaedle)
#     Benjamin M. Adams
#
# Further changes are done by:
#

############################
import numpy as np

def getIntersection(func1, func2, allowExtrapolation):

    #Both functions have to have only two columns
    if len(func1[0,:]) != 2:
        raise Exception('GetIntersection:OnlyTwoColumns - Function1 does not have two columns')

    if len(func2[0,:]) != 2:
        raise Exception('GetIntersection:OnlyTwoColumns - Function2 does not have two columns')

    #Both functions have to have at least two rows
    if len(func1[:,0]) < 2 or len(func2[:,0]) < 2:
        raise Exception('GetIntersection:NotEnoughRows - Functions do not both have more than two rows')

    # data points that directly cross
    data_int = []
    # data points where the extrapolation of one cross another
    extrap_data_int = []
    # data points where two extrapolations cross
    extrap_int = []

    # Iterate through all possible line combinations
    for index1 in range(len(func1[:,0])-1):
        for index2 in range(len(func2[:,0])-1):
            # Now x2(index2) is between x1(index1 -1) and x1(index1)
            x1a = func1[index1,0]
            y1a = func1[index1,1]
            x1b = func1[index1+1,0]
            y1b = func1[index1+1,1]
            x2a = func2[index2,0]
            y2a = func2[index2,1]
            x2b = func2[index2+1,0]
            y2b = func2[index2+1,1]

            # Solve for the equation of each line
            m1 = (y1b - y1a) / (x1b - x1a)
            b1 = y1a - m1*x1a
            m2 = (y2b - y2a) / (x2b - x2a)
            b2 = y2a - m2*x2a
            # Solve for x where the lines intersect
            x_int = (b2 - b1) / (m1 - m2)
            y_int = m1*x_int + b1
            # If x_int is between x1a and x1b, it crosses here
            if x_int >= x1a and x_int < x1b and x_int >= x2a and x_int < x2b:
                # Found an intersection!
                data_int.append(np.array([x_int, y_int]))
            elif allowExtrapolation == True:
                # See if extrapolated ends intersect anywhere
                if index2 == 0 and np.sign(x_int-x2a)==np.sign(x2a-x2b) and x_int >= x1a and x_int < x1b:
                    #It intersects this extrapolation!
                    extrap_data_int.append(np.array([x_int, y_int]))
                elif index2 == len(func2[:,0])-2 and np.sign(x_int-x2b)==np.sign(x2b-x2a) and x_int >= x1a and x_int < x1b:
                    #It intersects the extrapolation!
                    extrap_data_int.append(np.array([x_int, y_int]))

                if index1 == 0:
                    if np.sign(x_int-x1a)==np.sign(x1a-x1b) or len(func1[:,0])==2:
                        if x_int >= x2a and x_int < x2b:
                            #it intersects the existing line!
                            extrap_data_int.append(np.array([x_int, y_int]))
                        elif index2 == 0 and (np.sign(x_int-x2a)==np.sign(x2a-x2b) or len(func2[:,0])==2):
                            #It intersects this extrapolation!
                            extrap_int.append(np.array([x_int, y_int]))
                        elif index2 == len(func2[:,0])-2 and (np.sign(x_int-x2b)==np.sign(x2b-x2a) or len(func2[:,0])==2):
                            #It intersects the extrapolation!
                            extrap_int.append(np.array([x_int, y_int]))

                elif index1 == len(func1[:,0])-2:
                    if np.sign(x_int-x1b)==np.sign(x1b-x1a):
                        if x_int >= x2a and x_int < x2b:
                            #it intersects the existing line!
                            extrap_data_int.append(np.array([x_int, y_int]))
                        elif index2 == 0 and (np.sign(x_int-x2a)==np.sign(x2a-x2b) or len(func2[:,0])==2):
                            #It intersects this extrapolation!
                            extrap_int.append(np.array([x_int, y_int]))
                        elif index2 == len(func2[:,0])-2 and (np.sign(x_int-x2b)==np.sign(x2b-x2a) or len(func2[:,0])==2):
                            #It intersects the extrapolation!
                            extrap_int.append(np.array([x_int, y_int]))


    # Now we choose which of the results we want
    if data_int != []:
        # take first
        x_intercept = np.array(data_int)[0,0]
        y_intercept = np.array(data_int)[0,1]
    elif extrap_data_int != []:
        # take first
        x_intercept = np.array(extrap_data_int)[0,0]
        y_intercept = np.array(extrap_data_int)[0,1]
    elif extrap_int != []:
        # find closest value
        x_dx_min = np.inf
        extrap_int = np.array(extrap_int)
        for r in range(len(extrap_int[:,0])):
            x_int_left_dx = np.max([np.max(func1[:,0]), np.max(func2[:,0])]) - extrap_int[r,0]
            if x_int_left_dx < 0:
                x_int_left_dx = np.inf
            x_int_right_dx = extrap_int[r,0] - np.min([np.min(func1[:,0]), np.min(func2[:,0])])
            if x_int_right_dx < 0:
                x_int_right_dx = np.inf
            x_int_dx = np.min([x_int_left_dx, x_int_right_dx])
            if x_int_dx < x_dx_min:
                x_dx_min = x_int_dx
                # use this one
                x_intercept = extrap_int[r,0]
                y_intercept = extrap_int[r,1]

        #make sure x_dx_min was set
        if np.isinf(x_dx_min):
            raise Exception('GetIntersection:NoIntersection - Functions do not intersect')
    else:
        # Couldn't find any intersections
        raise Exception('GetIntersection:NoIntersection - Functions do not intersect')

    return (x_intercept, y_intercept)
