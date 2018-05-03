from tkinter import *
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
#from pprint import pprint

#TODO:
#   Bulletproof anti crash checking (move from bool about loading/supports to checking array)
#   Add shear/moment/deflection diagrams
#   Add editing of loading, supports, and beam length
#   distributed load across support

class Application(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.master = master

        self.SupportLocationArray = []
        self.SupportTypeArray = []
        self.BeamLength = 1
        self.BeamLength = 1
        self.PointLoadArray = []
        self.DistributedLoadArray = []
        self.IsBeamDefined = FALSE
        self.IsLoadDefined = FALSE
        self.IsSupportDefined = FALSE

        self.stillInsideLoad = FALSE
        self.xInsideLoad = 0

        plt.ion()

        self.ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=3)
        self.ax2 = plt.subplot2grid((3, 3), (2, 0), colspan=3)

        self.ax1.minorticks_on()
        self.ax1.grid(b=True, which='major', color='k', linestyle='-', alpha=0.3)
        self.ax1.grid(b=True, which='minor', color='g', linestyle='-', alpha=0.1)
        self.ax1.set_ylabel('V (kips)')
        self.ax1.set_xlabel('L (in)')

        self.ax2.minorticks_on()
        self.ax2.grid(b=True, which='major', color='k', linestyle='-', alpha=0.3)
        self.ax2.grid(b=True, which='minor', color='g', linestyle='-', alpha=0.1)
        self.ax2.set_ylabel('M (kips*in)')
        self.ax2.set_xlabel('L (in)')

        plt.show()

        self.init()

    def init(self):
        self.master.title("Beam Analysis")
        self.pack(fill=BOTH,expand=1)
        self.pack()
        menu = Menu(self.master)
        self.master.config(menu=menu)
        self.Canvas = Canvas(self.master, width=700, height=500)
        file = Menu(menu)

        file.add_command(label="Open")
        file.add_command(label="Save")
        file.add_command(label="Save As")
        file.add_command(label="About")
        file.add_command(label="Exit", command=self.Main_Exit)

        menu.add_cascade(label="File", menu=file)

        edit = Menu(menu)

        edit.add_command(label="Add Beam", command=self.AddBeam)
        edit.add_command(label="Add Support", command=self.AddSupport)
        edit.add_command(label="Add Loading", command=self.AddLoad)

        menu.add_cascade(label="Edit", menu=edit)

    def Main_Exit(self):
        exit()

    def validate(self, action, index, value_if_allowed, prior_value, text, validation_type, trigger_type, widget_name):
        if (action == '1'):
            if text in '0123456789.':
                try:
                    float(value_if_allowed)
                    return True
                except ValueError:
                    return False
            else:
                return False
        else:
            return True

    def AddBeam(self):
        self.AddBeam_Window = Toplevel(self)
        self.AddBeam_Window.geometry("250x200")
        self.AddBeamCanvas = Canvas(self.AddBeam_Window, width=200, height=200)

        Label(self.AddBeam_Window, text="Length (ft):").place(x=50, y=50, anchor="center")
        Label(self.AddBeam_Window, text="MOI (in^4):").place(x=50, y=100, anchor="center")

        vcmd = (self.master.register(self.validate), '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')

        self.BeamLengthEntry = Entry(self.AddBeam_Window, validate='key', validatecommand=vcmd, width=10)
        self.BeamMOIEntry = Entry(self.AddBeam_Window, validate='key', validatecommand=vcmd, width=10)

        self.BeamLengthEntry.place(x=175, y=50, anchor="center")
        self.BeamMOIEntry.place(x=175, y=100, anchor="center")

        Button(self.AddBeam_Window, text="OK", command=self.AddBeam_Exit, width=10).place(x=125, y=150, anchor="center")

        self.AddBeamCanvas.pack(fill=BOTH, expand=1)

    def AddSupport(self):

        if self.IsBeamDefined == TRUE:
            self.AddSupport_Window = Toplevel(self)
            self.AddSupport_Window.geometry("250x200")
            self.AddSupportCanvas = Canvas(self.AddSupport_Window, width=250, height=200)

            Label(self.AddSupport_Window, text="Location (ft):").place(x=50, y=50, anchor="center")

            vcmd = (self.master.register(self.validate), '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
            self.SupportLocationEntry = Entry(self.AddSupport_Window, validate='key', validatecommand=vcmd, width=10)

            self.SupportLocationEntry.place(x=175, y=50, anchor="center")

            self.Support = StringVar(self.master)
            self.SupportOptions = {"Fixed", "Roller", "Simple"}
            self.Support.set("Simple")
            OptionMenu(self.AddSupport_Window, self.Support, *self.SupportOptions).place(x=125, y=100, anchor="center")
            Button(self.AddSupport_Window, text="OK", command=self.AddSupport_Exit, width=10).place(x=125, y=150, anchor="center")

            self.AddSupportCanvas.pack(fill=BOTH, expand=1)

    def AddLoad(self):

        if self.IsBeamDefined == TRUE and self.IsSupportDefined == TRUE:
            self.AddLoad_Window = Toplevel(self)
            self.AddLoad_Window.geometry("150x100")

            Button(self.AddLoad_Window, text="Point Load", command=self.AddPointLoad).place(x=75, y=33, anchor="center")
            Button(self.AddLoad_Window, text="Distributed Load", command=self.AddDistributedLoad).place(x=75, y=66, anchor="center")

    def AddPointLoad(self):
        self.AddLoad_Window.destroy()
        self.AddPointLoad_Window = Toplevel(self)
        self.AddPointLoad_Window.geometry("275x200")
        self.AddPointLoadCanvas = Canvas(self.AddPointLoad_Window, width=275, height=200)

        Label(self.AddPointLoad_Window, text="Location (ft):").place(x=75, y=50, anchor="center")
        Label(self.AddPointLoad_Window, text="Magnitude (kip):").place(x=75, y=100, anchor="center")

        vcmd = (self.master.register(self.validate), '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')

        self.PointLoadLocationEntry = Entry(self.AddPointLoad_Window, validate='key', validatecommand=vcmd, width=10)
        self.PointLoadMagnitudeEntry = Entry(self.AddPointLoad_Window, validate='key', validatecommand=vcmd, width=10)

        self.PointLoadLocationEntry.place(x=200, y=50, anchor="center")
        self.PointLoadMagnitudeEntry.place(x=200, y=100, anchor="center")

        Button(self.AddPointLoad_Window, text="OK", command=self.AddPointLoad_Exit, width=10).place(x=137.5, y=150, anchor="center")

        self.AddPointLoadCanvas.pack(fill=BOTH, expand=1)

    def AddDistributedLoad(self):
        self.AddLoad_Window.destroy()
        self.AddDistributedLoad_Window = Toplevel(self)
        self.AddDistributedLoad_Window.geometry("325x325")
        self.AddDistributedLoadCanvas = Canvas(self.AddDistributedLoad_Window, width=300, height=325)

        Label(self.AddDistributedLoad_Window, text="Start location (ft):").place(x=100, y=50, anchor="center")
        Label(self.AddDistributedLoad_Window, text="Start magnitude (klf):").place(x=100, y=100, anchor="center")
        Label(self.AddDistributedLoad_Window, text="End location (ft):").place(x=100, y=175, anchor="center")
        Label(self.AddDistributedLoad_Window, text="End magnitude (klf):").place(x=100, y=225, anchor="center")

        vcmd = (self.master.register(self.validate), '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')

        self.DistributedLoadStartLocationEntry = Entry(self.AddDistributedLoad_Window, validate='key', validatecommand=vcmd, width=10)
        self.DistributedLoadStartMagnitudeEntry = Entry(self.AddDistributedLoad_Window, validate='key', validatecommand=vcmd, width=10)
        self.DistributedLoadEndLocationEntry = Entry(self.AddDistributedLoad_Window, validate='key', validatecommand=vcmd, width=10)
        self.DistributedLoadEndMagnitudeEntry = Entry(self.AddDistributedLoad_Window, validate='key', validatecommand=vcmd, width=10)

        self.DistributedLoadStartLocationEntry.place(x=225, y=50, anchor="center")
        self.DistributedLoadStartMagnitudeEntry.place(x=225, y=100, anchor="center")
        self.DistributedLoadEndLocationEntry.place(x=225, y=175, anchor="center")
        self.DistributedLoadEndMagnitudeEntry.place(x=225, y=225, anchor="center")

        Button(self.AddDistributedLoad_Window, text="OK", command=self.AddDistributedLoad_Exit, width=10).place(x=162.5, y=275, anchor="center")

        self.AddDistributedLoadCanvas.pack(fill=BOTH, expand=1)

    def AddBeam_Exit(self):

        if (float(self.BeamLengthEntry.get()) > 0) and (float(self.BeamMOIEntry.get()) > 0):
            self.BeamLength = float(self.BeamLengthEntry.get())
            self.BeamMOI = float(self.BeamMOIEntry.get())
            self.AddBeam_Window.destroy()
            self.Canvas.create_rectangle(50, 75, 650, 90, fill="#1f1", width=2)
            self.Canvas.create_line(50, 115, 50, 140, width=2)
            self.Canvas.create_line(650, 115, 650, 140, width=2)
            self.Canvas.create_line(325, 130, 50, 130, arrow=LAST, width=2)
            self.Canvas.create_line(375, 130, 650, 130, arrow=LAST, width=2)
            Label(self.master, text=str(self.BeamLength) + "'").place(x=350, y=130, anchor="center")
            self.IsBeamDefined = TRUE
            self.Canvas.pack(fill=BOTH, expand=1)

    def AddSupport_Exit(self):

        if (float(self.SupportLocationEntry.get()) >= 0) and (float(self.SupportLocationEntry.get()) <= float(self.BeamLength)):
            firstPoint = ((((float(self.SupportLocationEntry.get()))/float(self.BeamLength))*600)+50)

            if self.Support.get() == "Simple":
                self.SupportLocationArray.append(float(self.SupportLocationEntry.get()))
                self.SupportTypeArray.append("Simple")
                points = [firstPoint,90, firstPoint-10,105,firstPoint+10,105]
                self.Canvas.create_polygon(points, fill='red', width=1, outline='black')

            elif self.Support.get() == "Roller":
                self.SupportLocationArray.append(float(self.SupportLocationEntry.get()))
                self.SupportTypeArray.append("Roller")
                self.Canvas.create_oval(firstPoint-7.5,90,firstPoint+7.5,105, fill='red', width=1)

            elif self.Support.get() == "Fixed":
                if float(self.SupportLocationEntry.get()) == 0:
                    self.SupportLocationArray.append(float(self.SupportLocationEntry.get()))
                    self.SupportTypeArray.append("Fixed")

                    self.Canvas.create_line(50, 67.5, 50, 97.5, width=1)
                    self.Canvas.create_line(45, 72.5, 50, 67.5, width=1)
                    self.Canvas.create_line(45, 77.5, 50, 72.5, width=1)
                    self.Canvas.create_line(45, 82.5, 50, 77.5, width=1)
                    self.Canvas.create_line(45, 87.5, 50, 82.5, width=1)
                    self.Canvas.create_line(45, 92.5, 50, 87.5, width=1)
                    self.Canvas.create_line(45, 97.5, 50, 92.5, width=1)
                    self.Canvas.create_line(45, 102.5, 50, 97.5, width=1)

            self.Canvas.pack(fill=BOTH, expand=1)
            self.IsSupportDefined = TRUE
            self.AddSupport_Window.destroy()

    def AddPointLoad_Exit(self):

        if float(self.PointLoadMagnitudeEntry.get()) > 0:
            if (float(self.PointLoadLocationEntry.get()) > 0) and (float(self.PointLoadLocationEntry.get()) <= float(self.BeamLength)):

                entryPoint = ((((float(self.PointLoadLocationEntry.get())) / float(self.BeamLength)) * 600) + 50)

                newPointLoadArray = [float(self.PointLoadLocationEntry.get()), float(self.PointLoadMagnitudeEntry.get())]

                self.PointLoadArray.append(newPointLoadArray)

                self.Canvas.create_line(entryPoint, 35, entryPoint, 75, arrow=LAST, width=2, fill='blue')
                Label(self.master, text=self.PointLoadMagnitudeEntry.get() + " kip").place(x=entryPoint, y=25, anchor="center")

                self.Canvas.pack(fill=BOTH, expand=1)
                self.IsLoadDefined = TRUE
                self.AddPointLoad_Window.destroy()

    def AddDistributedLoad_Exit(self):

        if (float(self.DistributedLoadStartMagnitudeEntry.get()) > 0) and (float(self.DistributedLoadEndMagnitudeEntry.get()) > 0):
            if (float(self.DistributedLoadStartLocationEntry.get()) >= 0) and (float(self.DistributedLoadEndLocationEntry.get()) <= float(self.BeamLength)):
                if (float(self.DistributedLoadEndLocationEntry.get()) >= 0) and (float(self.DistributedLoadStartLocationEntry.get()) <= float(self.BeamLength)):
                    if (float(self.DistributedLoadEndLocationEntry.get())) > (float(self.DistributedLoadStartLocationEntry.get())):

                        startEntryPoint = ((((float(self.DistributedLoadStartLocationEntry.get())) / float(self.BeamLength)) * 600) + 50)
                        endEntryPoint = ((((float(self.DistributedLoadEndLocationEntry.get())) / float(self.BeamLength)) * 600) + 50)

                        newDistributedLoadArray = [float(self.DistributedLoadStartLocationEntry.get()), float(self.DistributedLoadStartMagnitudeEntry.get()), float(self.DistributedLoadEndLocationEntry.get()), float(self.DistributedLoadEndMagnitudeEntry.get())]

                        self.DistributedLoadArray.append(newDistributedLoadArray)

                        x=0
                        while (startEntryPoint+(15*x)) <= endEntryPoint:
                            self.Canvas.create_line(startEntryPoint+(15*x), 35, startEntryPoint+(15*x), 75, arrow=LAST, width=2, fill='blue')
                            x=x+1

                        self.Canvas.create_line(startEntryPoint, 35, endEntryPoint, 35, width=2, fill='blue')
                        Label(self.master, text=self.DistributedLoadStartMagnitudeEntry.get() + " klf").place(x=startEntryPoint, y=20, anchor="center")
                        Label(self.master, text=self.DistributedLoadEndMagnitudeEntry.get() + " klf").place(x=endEntryPoint, y=20, anchor="center")

                        self.Canvas.pack(fill=BOTH, expand=1)
                        self.IsLoadDefined = TRUE
                        self.AddDistributedLoad_Window.destroy()

    def Analysis(self):

        #Stability checks
        if (self.IsBeamDefined == TRUE) and (self.IsLoadDefined == TRUE) and (self.IsSupportDefined == TRUE):

            #Setup element variables
            spanFixedEndMatrix = []
            spanForceMatrix = []
            spanStiffnessMatrix = []

            #Align support location array
            sortedSupportArray = sorted(zip(self.SupportLocationArray, self.SupportTypeArray))
            self.SupportLocationArray = [item[0] for item in sortedSupportArray]
            self.SupportTypeArray = [item[1] for item in sortedSupportArray]

            sortedPointLoadArray = sorted(self.PointLoadArray)
            self.PointLoadArray = sortedPointLoadArray

            sortedDistributedLoadArray = sorted(self.DistributedLoadArray)
            self.DistributedLoadArray = sortedDistributedLoadArray

            #Add to support location array if any cantilever spans are encountered
            if (float(self.SupportLocationArray[0]) != 0.0):
                self.SupportLocationArray.insert(0, 0.0)
                self.SupportTypeArray.insert(0, "Free")
            elif (float(self.SupportLocationArray[len(self.SupportLocationArray) - 1]) != float(self.BeamLength)):
                self.SupportLocationArray.append(self.BeamLength)
                self.SupportTypeArray.append("Free")

            #More stability checks
            if (len(self.SupportLocationArray) >= 2) or (self.SupportTypeArray[0] == "Fixed"):

                #Calculate number of spans to loop through
                if (len(self.SupportLocationArray) >= 2):
                    numberSpans = (len(self.SupportLocationArray)) - 1
                else:
                    numberSpans = 1

                for x in range(0, numberSpans):

                    #Calculate each span length
                    if numberSpans > 1:
                        spanLength = float(self.SupportLocationArray[x+1]) - float(self.SupportLocationArray[x])
                    else:
                        spanLength = float(self.BeamLength)

                    #Initialize element variables
                    spanFixedEndMatrix.append(np.array([[0.0], [0.0], [0.0], [0.0]]))
                    spanForceMatrix.append(np.array([[0.0], [0.0], [0.0], [0.0]]))
                    spanStiffnessMatrix.append(np.array([[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]]))

                    #Loop through point load array
                    for y in range(0, len(self.PointLoadArray)):
                        pointLoadLocationArray = [item[0] for item in self.PointLoadArray]
                        pointLoadMagnitudeArray = [item[1] for item in self.PointLoadArray]

                        #More stability checks
                        if (pointLoadLocationArray[y] <= float(self.SupportLocationArray[x+1])) and (pointLoadLocationArray[y] >= float(self.SupportLocationArray[x])):

                            #Calculate fixed end forces for any point loads on current span and add them to element fixed end force matrix
                            a=(float(self.SupportLocationArray[x+1])-pointLoadLocationArray[y])*12
                            b=((spanLength*12)-a)

                            m1 = float((pointLoadMagnitudeArray[y]*(b**2)*a)/((spanLength*12)**2))
                            m2 = float((pointLoadMagnitudeArray[y]*(a**2)*b)/((spanLength*12)**2))
                            r1 = float(((pointLoadMagnitudeArray[y]*(b**2))/((spanLength*12)**3))*((3*a)+b))
                            r2 = float((pointLoadMagnitudeArray[y]*(a**2)*(a+(3*b)))/((spanLength*12)**3))

                            spanFixedEndMatrix[x][0] = spanFixedEndMatrix[x][0] - r2
                            spanFixedEndMatrix[x][1] = spanFixedEndMatrix[x][1] - m2
                            spanFixedEndMatrix[x][2] = spanFixedEndMatrix[x][2] - r1
                            spanFixedEndMatrix[x][3] = spanFixedEndMatrix[x][3] + m1

                    #Loop through distributed load array
                    for z in range(0, len(self.DistributedLoadArray)):
                        distributedLoadStartLocationArray = [item[0] for item in self.DistributedLoadArray]
                        distributedLoadEndLocationArray = [item[2] for item in self.DistributedLoadArray]
                        distributedLoadStartMagnitudeArray = [item[1] for item in self.DistributedLoadArray]
                        distributedLoadEndMagnitudeArray = [item[3] for item in self.DistributedLoadArray]

                        #More stability checks
                        if (distributedLoadEndLocationArray[z] <= float(self.SupportLocationArray[x+1])) and (distributedLoadStartLocationArray[z] >= float(self.SupportLocationArray[x])):

                            #If uniform partial or fully distibuted load, calculate fixed end forces and add them to element fixed end force matrix
                            if (distributedLoadStartMagnitudeArray[z]==distributedLoadEndMagnitudeArray[z]):

                                a1=(distributedLoadStartLocationArray[z]-float(self.SupportLocationArray[x]))*12
                                b1=((distributedLoadEndLocationArray[z]-distributedLoadStartLocationArray[z])*12)+a1
                                c1=((distributedLoadEndLocationArray[z]-distributedLoadStartLocationArray[z])*12)
                                d1=(spanLength*12)-(a1/2)-(b1/2)
                                L=spanLength*12

                                m11 = float((((-distributedLoadStartMagnitudeArray[z]/12)*c1)/(24*L))*(((24*(d1**3))/L) - ((6*b1*(c1**2))/L) + ((3*(c1**3))/L) + (4*(c1**2)) - (24*(d1**2))))
                                m21 = float((((distributedLoadStartMagnitudeArray[z]/12)*c1)/(24*L))*(((24*(d1**3))/L) - ((6*b1*(c1**2))/L) + ((3*(c1**3))/L) + (2*(c1**2)) - (48*(d1**2)) + (24*d1*L)))
                                r11 = float((((distributedLoadStartMagnitudeArray[z]/12)*c1)/(4*(L**2)))*((12*(d1**2)) - ((8*(d1**3))/L) + ((2*b1*(c1**2))/L) - ((c1**3)/L) - (c1**2)))
                                r21 = float(((distributedLoadStartMagnitudeArray[z]/12)*c1)-r11)

                                spanFixedEndMatrix[x][0] = spanFixedEndMatrix[x][0] - r11
                                spanFixedEndMatrix[x][1] = spanFixedEndMatrix[x][1] - m11
                                spanFixedEndMatrix[x][2] = spanFixedEndMatrix[x][2] - r21
                                spanFixedEndMatrix[x][3] = spanFixedEndMatrix[x][3] + m21

                            #If increasing partial or fully distributed load
                            elif (distributedLoadStartMagnitudeArray[z] < distributedLoadEndMagnitudeArray[z]):
                                a1 = (distributedLoadStartLocationArray[z] - float(self.SupportLocationArray[x])) * 12
                                b1 = ((distributedLoadEndLocationArray[z] - distributedLoadStartLocationArray[z]) * 12) + a1
                                c1 = ((distributedLoadEndLocationArray[z] - distributedLoadStartLocationArray[z]) * 12)
                                d1 = (spanLength * 12) - (a1 / 2) - (b1 / 2)
                                d2 = (spanLength * 12) - a1 - ((2*c1)/3)
                                w = (distributedLoadEndMagnitudeArray[z] - distributedLoadStartMagnitudeArray[z])/12
                                W = (w*(c1))/2
                                L = spanLength * 12

                                m11 = float((((-distributedLoadStartMagnitudeArray[z]/12) * c1) / (24 * L)) * (((24 * (d1 ** 3)) / L) - ((6 * b1 * (c1 ** 2)) / L) + ((3 * (c1 ** 3)) / L) + (4*(c1 ** 2)) - (24 * (d1 ** 2))))
                                m21 = float((((distributedLoadStartMagnitudeArray[z]/12) * c1) / (24 * L)) * (((24 * (d1 ** 3)) / L) - ((6 * b1 * (c1 ** 2)) / L) + ((3 * (c1 ** 3)) / L) + (2 * (c1 ** 2)) - (48 * (d1 ** 2)) + (24 * d1 * L)))
                                r11 = float((((distributedLoadStartMagnitudeArray[z]/12) * c1) / (4 * (L ** 2))) * ((12 * (d1 ** 2)) - ((8 * (d1 ** 3)) / L) + ((2 * b1 * (c1 ** 2)) / L) - ((c1 ** 3) / L) - (c1 ** 2)))
                                r21 = float(((distributedLoadStartMagnitudeArray[z]/12) * c1) - r11)

                                m12 = float(((-w*c1)/(2*(L**2)))*( ((d2**2)*(d2-L)) + (((c1**2)/3)*((L/3)+((17/90)*c1)-(b1/2)))))
                                m22 = float(((w * c1) / (2 * (L ** 2))) * (((d2)*((d2 - L)**2)) + (((c1 ** 2) / 6)*((L / 3) + ((17 / 45) * c1) - (b1)))))
                                r12 = float(((w*c1)/(2*(L**3)))*( ((d2**2)*((3*L)-(2*d2))) - (((c1**2)/3)*((L/2)-b1+((17/45)*c1)))))
                                r22 = float(W-r12)

                                spanFixedEndMatrix[x][0] = spanFixedEndMatrix[x][0] - (r11 + r12)
                                spanFixedEndMatrix[x][1] = spanFixedEndMatrix[x][1] - (m11 + m12)
                                spanFixedEndMatrix[x][2] = spanFixedEndMatrix[x][2] - (r21 + r22)
                                spanFixedEndMatrix[x][3] = spanFixedEndMatrix[x][3] + (m21 + m22)

                            #If decreasing partial or fully distibuted load
                            elif (distributedLoadStartMagnitudeArray[z] > distributedLoadEndMagnitudeArray[z]):
                                a1 = (distributedLoadStartLocationArray[z] - float(self.SupportLocationArray[x])) * 12
                                a2 = (float(self.SupportLocationArray[x+1]) - distributedLoadEndLocationArray[z]) * 12
                                b1 = ((distributedLoadEndLocationArray[z] - distributedLoadStartLocationArray[z]) * 12) + a1
                                b2 = ((distributedLoadEndLocationArray[z] - distributedLoadStartLocationArray[z]) * 12) + a2
                                c1 = ((distributedLoadEndLocationArray[z] - distributedLoadStartLocationArray[z]) * 12)
                                d1 = (spanLength * 12) - (a1 / 2) - (b1 / 2)
                                d2 = (spanLength * 12) - a2 - ((2 * c1) / 3)
                                w = (distributedLoadStartMagnitudeArray[z] - distributedLoadEndMagnitudeArray[z])/12
                                W = (w * (c1)) / 2
                                L = spanLength * 12

                                m11 = float((((-distributedLoadEndMagnitudeArray[z]/12) * c1) / (24 * L)) * (((24 * (d1 ** 3)) / L) - ((6 * b1 * (c1 ** 2)) / L) + ((3 * (c1 ** 3)) / L) + (4*(c1 ** 2)) - (24 * (d1 ** 2))))
                                m21 = float((((distributedLoadEndMagnitudeArray[z]/12) * c1) / (24 * L)) * (((24 * (d1 ** 3)) / L) - ((6 * b1 * (c1 ** 2)) / L) + ((3 * (c1 ** 3)) / L) + (2 * (c1 ** 2)) - (48 * (d1 ** 2)) + (24 * d1 * L)))
                                r11 = float((((distributedLoadEndMagnitudeArray[z]/12) * c1) / (4 * (L ** 2))) * ((12 * (d1 ** 2)) - ((8 * (d1 ** 3)) / L) + ((2 * b1 * (c1 ** 2)) / L) - ((c1 ** 3) / L) - (c1 ** 2)))
                                r21 = float(((distributedLoadEndMagnitudeArray[z]/12) * c1) - r11)

                                m12 = float(((-w * c1) / (2 * (L ** 2))) * (((d2 ** 2)*(d2 - L)) + (((c1 ** 2) / 3)*((L / 3) + ((17 / 90) * c1) - (b2 / 2)))))
                                m22 = float(((w * c1) / (2 * (L ** 2))) * (((d2)*((d2 - L) ** 2)) + (((c1 ** 2) / 6)*((L / 3) + ((17 / 45) * c1) - (b2)))))
                                r12 = float(((w * c1) / (2 * (L ** 3))) * (((d2 ** 2)*((3 * L) - (2 * d2))) - (((c1 ** 2) / 3)*((L / 2) - b2 + ((17 / 45) * c1)))))
                                r22 = float(W - r12)

                                spanFixedEndMatrix[x][0] = spanFixedEndMatrix[x][0] - (r11 + r22)
                                spanFixedEndMatrix[x][1] = spanFixedEndMatrix[x][1] - (m11 + m22)
                                spanFixedEndMatrix[x][2] = spanFixedEndMatrix[x][2] - (r21 + r12)
                                spanFixedEndMatrix[x][3] = spanFixedEndMatrix[x][3] + (m21 + m12)

                    E = 29000
                    I = float(self.BeamMOI)
                    L = spanLength*12

                    spanStiffnessMatrix[x]=np.array([[((E*I)/(L**3))*12, ((E*I)/(L**3))*6*L, ((E*I)/(L**3))*-12, ((E*I)/(L**3))*6*L ],
                                             [((E*I)/(L**3))*6*L,((E*I)/(L**3))*4*(L**2), ((E*I)/(L**3))*-6*L,((E*I)/(L**3))*2*(L**2)],
                                             [((E*I)/(L**3))*-12,((E*I)/(L**3))*-6*L, ((E*I)/(L**3))*12,((E*I)/(L**3))*-6*L],
                                             [((E*I)/(L**3))*6*L,((E*I)/(L**3))*2*(L**2),((E*I)/(L**3))*-6*L,((E*I)/(L**3))*4*(L**2)]])


            if len(spanForceMatrix) > 1:
                globalStiffnessMatrix = spanStiffnessMatrix[0]
                globalFixedEndMatrix = spanFixedEndMatrix[0]
                dispMatrix = []
                globalForceMatrix = []

                for x in range(1, len(spanForceMatrix)):
                    col1 = []
                    col2 = []
                    row1 = []
                    row2 = []

                    globalStiffnessMatrix[(2*x), (2*x)] = globalStiffnessMatrix[(2*x), (2*x)] + spanStiffnessMatrix[x][0,0]
                    globalStiffnessMatrix[(2*x)+1, (2*x)] = globalStiffnessMatrix[(2*x)+1, (2*x)] + spanStiffnessMatrix[x][1,0]
                    globalStiffnessMatrix[(2*x), (2*x)+1] = globalStiffnessMatrix[(2*x), (2*x)+1] + spanStiffnessMatrix[x][0,1]
                    globalStiffnessMatrix[(2*x)+1, (2*x)+1] = globalStiffnessMatrix[(2*x)+1, (2*x)+1] + spanStiffnessMatrix[x][1,1]

                    for y in range(1, x+1):
                        col1.append(np.array(0.0))
                        col1.append(np.array(0.0))
                        col2.append(np.array(0.0))
                        col2.append(np.array(0.0))
                        row1.append(np.array(0.0))
                        row1.append(np.array(0.0))
                        row2.append(np.array(0.0))
                        row2.append(np.array(0.0))

                    for i in range(0,2):
                        col1.append(np.array(spanStiffnessMatrix[x][i, 2]))
                        col2.append(np.array(spanStiffnessMatrix[x][i, 3]))

                    for j in range(0,4):
                        row1.append(np.array(spanStiffnessMatrix[x][2,j]))
                        row2.append(np.array(spanStiffnessMatrix[x][3,j]))

                    globalStiffnessMatrix = np.column_stack([globalStiffnessMatrix, col1])
                    globalStiffnessMatrix = np.column_stack([globalStiffnessMatrix, col2])
                    globalStiffnessMatrix = np.row_stack([globalStiffnessMatrix, row1])
                    globalStiffnessMatrix = np.row_stack([globalStiffnessMatrix, row2])

                    row = []

                    globalFixedEndMatrix[(2*x)] = globalFixedEndMatrix[(2*x)] + spanFixedEndMatrix[x][0]
                    globalFixedEndMatrix[(2*x)+1] = globalFixedEndMatrix[(2*x)+1] + spanFixedEndMatrix[x][1]

                    row.append(np.array(spanFixedEndMatrix[x][2]))
                    row.append(np.array(spanFixedEndMatrix[x][3]))

                    globalFixedEndMatrix = np.row_stack([globalFixedEndMatrix, row])

                globalModifiedStiffnessMatrix = globalStiffnessMatrix
                globalModifiedFixedEndMatrix = globalFixedEndMatrix

                if (self.SupportTypeArray[0] == "Free"):

                    if (self.SupportTypeArray[1] == "Simple") or (self.SupportTypeArray[1] == "Roller"):
                        globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (2), axis=0)
                        globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (2), axis=1)
                        globalModifiedFixedEndMatrix = np.delete(globalModifiedFixedEndMatrix, (2), axis=0)

                    for x in range(1, len(spanForceMatrix)):

                        if (self.SupportTypeArray[x + 1] == "Simple") or (self.SupportTypeArray[x + 1] == "Roller"):
                            globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (x + 2), axis=0)
                            globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (x + 2), axis=1)
                            globalModifiedFixedEndMatrix = np.delete(globalModifiedFixedEndMatrix, (x + 2), axis=0)

                elif (self.SupportTypeArray[0] == "Fixed"):

                    globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (0), axis=0)
                    globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (0), axis=1)
                    globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (0), axis=0)
                    globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (0), axis=1)
                    globalModifiedFixedEndMatrix = np.delete(globalModifiedFixedEndMatrix, (0), axis=0)
                    globalModifiedFixedEndMatrix = np.delete(globalModifiedFixedEndMatrix, (0), axis=0)

                    if (self.SupportTypeArray[1] == "Simple") or (self.SupportTypeArray[1] == "Roller"):
                        globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (0), axis=0)
                        globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (0), axis=1)
                        globalModifiedFixedEndMatrix = np.delete(globalModifiedFixedEndMatrix, (0), axis=0)

                    for x in range(1, len(spanForceMatrix)):

                        if (self.SupportTypeArray[x + 1] == "Simple") or (self.SupportTypeArray[x + 1] == "Roller"):
                            globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (x), axis=0)
                            globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (x), axis=1)
                            globalModifiedFixedEndMatrix = np.delete(globalModifiedFixedEndMatrix, (x), axis=0)

                else:

                    globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (0), axis=0)
                    globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (0), axis=1)
                    globalModifiedFixedEndMatrix = np.delete(globalModifiedFixedEndMatrix, (0), axis=0)

                    if (self.SupportTypeArray[1] == "Simple") or (self.SupportTypeArray[1] == "Roller"):
                        globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (1), axis=0)
                        globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (1), axis=1)
                        globalModifiedFixedEndMatrix = np.delete(globalModifiedFixedEndMatrix, (1), axis=0)

                    for x in range(1, len(spanForceMatrix)):

                        if (self.SupportTypeArray[x + 1] == "Simple") or (self.SupportTypeArray[x + 1] == "Roller"):
                            globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (x+1), axis=0)
                            globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (x+1), axis=1)
                            globalModifiedFixedEndMatrix = np.delete(globalModifiedFixedEndMatrix, (x+1), axis=0)

            elif len(spanForceMatrix) == 1:

                    globalStiffnessMatrix = spanStiffnessMatrix[0]
                    globalFixedEndMatrix = spanFixedEndMatrix[0]
                    dispMatrix = []
                    globalForceMatrix = []

                    globalModifiedStiffnessMatrix = globalStiffnessMatrix
                    globalModifiedFixedEndMatrix = globalFixedEndMatrix

                    if (self.SupportTypeArray[0] == "Free"):

                        if (self.SupportTypeArray[1] == "Simple") or (self.SupportTypeArray[1] == "Roller"):
                            globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (2), axis=0)
                            globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (2), axis=1)
                            globalModifiedFixedEndMatrix = np.delete(globalModifiedFixedEndMatrix, (2), axis=0)

                    elif (self.SupportTypeArray[0] == "Fixed"):

                        globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (0), axis=0)
                        globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (0), axis=1)
                        globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (0), axis=0)
                        globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (0), axis=1)
                        globalModifiedFixedEndMatrix = np.delete(globalModifiedFixedEndMatrix, (0), axis=0)
                        globalModifiedFixedEndMatrix = np.delete(globalModifiedFixedEndMatrix, (0), axis=0)

                        if (self.SupportTypeArray[1] == "Simple") or (self.SupportTypeArray[1] == "Roller"):
                            globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (0), axis=0)
                            globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (0), axis=1)
                            globalModifiedFixedEndMatrix = np.delete(globalModifiedFixedEndMatrix, (0), axis=0)

                    else:

                        globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (0), axis=0)
                        globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (0), axis=1)
                        globalModifiedFixedEndMatrix = np.delete(globalModifiedFixedEndMatrix, (0), axis=0)

                        if (self.SupportTypeArray[1] == "Simple") or (self.SupportTypeArray[1] == "Roller"):
                            globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (1), axis=0)
                            globalModifiedStiffnessMatrix = np.delete(globalModifiedStiffnessMatrix, (1), axis=1)
                            globalModifiedFixedEndMatrix = np.delete(globalModifiedFixedEndMatrix, (1), axis=0)

            invGlobalStiffnessMatrix = np.linalg.inv(globalModifiedStiffnessMatrix)
            dispMatrix = invGlobalStiffnessMatrix.dot(np.negative(globalModifiedFixedEndMatrix))

            if len(spanForceMatrix) > 1:

                if (self.SupportTypeArray[0] == "Free"):

                    if (self.SupportTypeArray[1] == "Simple") or (self.SupportTypeArray[1] == "Roller"):
                        dispMatrix = np.insert(dispMatrix, (2), np.array(0.0), 0)

                    for x in range(1, numberSpans):

                        if (self.SupportTypeArray[x + 1] == "Simple") or (self.SupportTypeArray[x + 1] == "Roller"):
                            dispMatrix = np.insert(dispMatrix, ((2*x)+2), np.array(0.0), 0)

                elif (self.SupportTypeArray[0] == "Fixed"):

                    dispMatrix = np.insert(dispMatrix, (0), np.array(0.0), 0)
                    dispMatrix = np.insert(dispMatrix, (0), np.array(0.0), 0)

                    if (self.SupportTypeArray[1] == "Simple") or (self.SupportTypeArray[1] == "Roller"):
                        dispMatrix = np.insert(dispMatrix, (2), np.array(0.0), 0)

                    for x in range(1, numberSpans):

                        if (self.SupportTypeArray[x + 1] == "Simple") or (self.SupportTypeArray[x + 1] == "Roller"):
                            dispMatrix = np.insert(dispMatrix, ((2 * x)+2), np.array(0.0), 0)
                else:
                    dispMatrix = np.insert(dispMatrix, (0), np.array(0.0), 0)

                    if (self.SupportTypeArray[1] == "Simple") or (self.SupportTypeArray[1] == "Roller"):
                        dispMatrix = np.insert(dispMatrix, (2), np.array(0.0), 0)

                    for x in range(1, numberSpans):

                        if (self.SupportTypeArray[x + 1] == "Simple") or (self.SupportTypeArray[x + 1] == "Roller"):
                            dispMatrix = np.insert(dispMatrix, ((2 * x) + 2), np.array(0.0), 0)

            elif len(spanForceMatrix) == 1:
                if (self.SupportTypeArray[0] == "Simple") or (self.SupportTypeArray[0] == "Roller"):
                    dispMatrix = np.insert(dispMatrix, (0), np.array(0.0), 0)
                elif (self.SupportTypeArray[0] == "Fixed"):
                    dispMatrix = np.insert(dispMatrix, (0), np.array(0.0), 0)
                    dispMatrix = np.insert(dispMatrix, (0), np.array(0.0), 0)

                if (self.SupportTypeArray[1] == "Simple") or (self.SupportTypeArray[1] == "Roller"):
                    dispMatrix = np.insert(dispMatrix, (2), np.array(0.0), 0)

            globalForceMatrix = np.add(globalStiffnessMatrix.dot(dispMatrix), globalFixedEndMatrix)

            shearX = []
            shearY = []
            momentX = []
            momentY = []
            spanLength= []
            cumX = 0
            pointLoadIndexMoment = []
            distLoadIndexMoment = []

            for i in range(0, len(spanForceMatrix)):
                pointLoadIndex = []
                distLoadIndex = []

                if len(spanForceMatrix) > 1:
                    spanLength.append(float(self.SupportLocationArray[i + 1]) - float(self.SupportLocationArray[i]))

                else:
                    spanLength.append(float(self.BeamLength))

                spanLength[i]=spanLength[i]*12

                if (len(shearX))==0:
                    shearY.append(0.0)
                    shearX.append(0.0)

                if len(shearX) > 0:
                    prevSpanShear = shearY[int(cumX)]
                else:
                    prevSpanShear = 0.0

                for x in range(0, int(spanLength[i])+1):

                    pointLoadLocationArray = [item[0] for item in self.PointLoadArray]
                    pointLoadMagnitudeArray = [item[1] for item in self.PointLoadArray]

                    distributedLoadStartLocationArray = [item[0] for item in self.DistributedLoadArray]
                    distributedLoadEndLocationArray = [item[2] for item in self.DistributedLoadArray]
                    distributedLoadStartMagnitudeArray = [item[1] for item in self.DistributedLoadArray]
                    distributedLoadEndMagnitudeArray = [item[3] for item in self.DistributedLoadArray]

                    shearY[int(cumX+x)] = shearY[int(cumX+x)] + prevSpanShear + float(-1*globalForceMatrix[(2*i)])

                    if (float((cumX+x)/12) in pointLoadLocationArray) == TRUE:
                        index = pointLoadLocationArray.index(float((cumX+x)/12))
                        pointLoadIndex.append(index)

                        if index not in pointLoadIndexMoment:
                            pointLoadIndexMoment.append(index)

                    if ((float((cumX + x) / 12) in distributedLoadStartLocationArray) == TRUE):
                        index = distributedLoadStartLocationArray.index((float((cumX + x) / 12)))
                        distLoadIndex.append(index)

                        if index not in distLoadIndexMoment:
                            distLoadIndexMoment.append(index)

                    if len(pointLoadIndex) > 0:

                        for j in range(0, len(pointLoadIndex)):

                            if x > (x-((pointLoadLocationArray[pointLoadIndex[j]]*12)-cumX)):
                                shearY[int(cumX+x)] = shearY[int(cumX + x)] + (float(-1*pointLoadMagnitudeArray[pointLoadIndex[j]]*((x-((pointLoadLocationArray[pointLoadIndex[j]]*12)-cumX))**0)))

                    if len(distLoadIndex) > 0:

                        for j in range(0, len(distLoadIndex)):

                            loadLength = float(distributedLoadEndLocationArray[distLoadIndex[j]] - distributedLoadStartLocationArray[distLoadIndex[j]])*12.0

                            if (distributedLoadStartMagnitudeArray[distLoadIndex[j]] == distributedLoadEndMagnitudeArray[distLoadIndex[j]]):

                                if x > ((((distributedLoadStartLocationArray[distLoadIndex[j]]*12)-cumX))):
                                    shearY[int(cumX + x)] = shearY[int(cumX + x)]+(float((-1*distributedLoadStartMagnitudeArray[distLoadIndex[j]]/12)*((x-((distributedLoadStartLocationArray[distLoadIndex[j]]*12)-cumX))**1)))

                                if ((distributedLoadStartLocationArray[distLoadIndex[j]]*12)+loadLength-cumX) < spanLength[i]:
                                    if x > ((distributedLoadStartLocationArray[distLoadIndex[j]]*12)+loadLength-cumX):
                                        compX = ((distributedLoadStartLocationArray[distLoadIndex[j]]*12)+loadLength-cumX)
                                        shearY[int(cumX + x)] = shearY[int(cumX + x)] + (float(( distributedLoadStartMagnitudeArray[distLoadIndex[j]] / 12) * ((x - compX) ** 1)))

                            elif (distributedLoadStartMagnitudeArray[distLoadIndex[j]] < distributedLoadEndMagnitudeArray[distLoadIndex[j]]):

                                if ((distributedLoadStartLocationArray[distLoadIndex[j]]*12)+loadLength-cumX) < spanLength[i]:

                                    compX = ((distributedLoadStartLocationArray[distLoadIndex[j]]*12)+loadLength-cumX)
                                    w = (float((distributedLoadEndMagnitudeArray[distLoadIndex[j]]-distributedLoadStartMagnitudeArray[distLoadIndex[j]]))/12.0)
                                    k = (w/loadLength)

                                    if x > ((((distributedLoadStartLocationArray[distLoadIndex[j]] * 12) - cumX))):
                                        shearY[int(cumX + x)] = shearY[int(cumX + x)] + (float(((-1*distributedLoadStartMagnitudeArray[distLoadIndex[j]]/12) * ((x-((distributedLoadStartLocationArray[distLoadIndex[j]]*12)-cumX))**1))+(-1*(k/2)*((x-((distributedLoadStartLocationArray[distLoadIndex[j]]*12)-cumX))**2))))
                                    if x > ((distributedLoadStartLocationArray[distLoadIndex[j]] * 12) + loadLength - cumX):
                                        shearY[int(cumX + x)] = shearY[int(cumX + x)] + (float(((distributedLoadStartMagnitudeArray[distLoadIndex[j]]/12) * ((x-compX)**1))+((k/2)*((x-compX)**2)))) + (w*((x-compX)**1))
                                else:
                                    if x > ((((distributedLoadStartLocationArray[distLoadIndex[j]] * 12) - cumX))):
                                        k = float((distributedLoadEndMagnitudeArray[distLoadIndex[j]]-distributedLoadStartMagnitudeArray[distLoadIndex[j]])/loadLength)/12.0
                                        shearY[int(cumX + x)] = shearY[int(cumX + x)] + (float(((-1*distributedLoadStartMagnitudeArray[distLoadIndex[j]]/12) * ((x-((distributedLoadStartLocationArray[distLoadIndex[j]]*12)-cumX))**1))+(-1*(k/2)*((x-((distributedLoadStartLocationArray[distLoadIndex[j]]*12)-cumX))**2))))

                            elif (distributedLoadStartMagnitudeArray[distLoadIndex[j]] > distributedLoadEndMagnitudeArray[distLoadIndex[j]]):

                                if ((distributedLoadStartLocationArray[distLoadIndex[j]]*12)+loadLength-cumX) < spanLength[i]:

                                    compX = ((distributedLoadStartLocationArray[distLoadIndex[j]]*12)+loadLength-cumX)
                                    w = (float((distributedLoadStartMagnitudeArray[distLoadIndex[j]]-distributedLoadEndMagnitudeArray[distLoadIndex[j]]))/12.0)
                                    k = (w/loadLength)

                                    if x > ((((distributedLoadStartLocationArray[distLoadIndex[j]] * 12) - cumX))):
                                        shearY[int(cumX + x)] = shearY[int(cumX + x)]+ (float(((-1*distributedLoadEndMagnitudeArray[distLoadIndex[j]]/12) * ((x-((distributedLoadStartLocationArray[distLoadIndex[j]]*12)-cumX))**1))-(-1*(k/2)*((x-((distributedLoadStartLocationArray[distLoadIndex[j]]*12)-cumX))**2))+((-1*w) * ((x-((distributedLoadStartLocationArray[distLoadIndex[j]]*12)-cumX))**1))))
                                    if x > ((distributedLoadStartLocationArray[distLoadIndex[j]] * 12) + loadLength - cumX):
                                        shearY[int(cumX + x)] = shearY[int(cumX + x)] + (float(((distributedLoadEndMagnitudeArray[distLoadIndex[j]]/12) * ((x-compX)**1))-((k/2)*(((x)-compX)**2)) - (w*(((x)-compX)**1))+((w) * ((x-compX)**1))))
                                else:
                                    if x > ((((distributedLoadStartLocationArray[distLoadIndex[j]] * 12) - cumX))):
                                        k = float((distributedLoadStartMagnitudeArray[distLoadIndex[j]]-distributedLoadEndMagnitudeArray[distLoadIndex[j]])/loadLength)/12.0
                                        shearY[int(cumX + x)] = shearY[int(cumX + x)] + (float(((-1*distributedLoadEndMagnitudeArray[distLoadIndex[j]]/12) * ((x-((distributedLoadStartLocationArray[distLoadIndex[j]]*12)-cumX))**1))+(-1*(k/2)*((x-((distributedLoadStartLocationArray[distLoadIndex[j]]*12)-cumX))**2))))

                    if x != spanLength[i]:
                        shearX.append(float(cumX + x + 1))
                        shearY.append(0.0)

                if (len(momentX))==0:
                    if (self.SupportTypeArray[0] != "Free"):
                        momentY.append(-1*globalForceMatrix[(2 * i) + 1])
                    else:
                        momentY.append(0.0)
                    momentX.append(0.0)

                for x in range(0, int(spanLength[i])+1):

                    pointLoadLocationArray = [item[0] for item in self.PointLoadArray]
                    pointLoadMagnitudeArray = [item[1] for item in self.PointLoadArray]

                    distributedLoadStartLocationArray = [item[0] for item in self.DistributedLoadArray]
                    distributedLoadEndLocationArray = [item[2] for item in self.DistributedLoadArray]
                    distributedLoadStartMagnitudeArray = [item[1] for item in self.DistributedLoadArray]
                    distributedLoadEndMagnitudeArray = [item[3] for item in self.DistributedLoadArray]

                    if x > 0:
                        momentY[int(cumX + x)] = momentY[int(cumX + x)] + float(-1*(globalForceMatrix[(2 * i) + 1])) + float((float(globalForceMatrix[(2 * i)]) * x))

                        for k in range (0,i):
                            momentY[int(cumX + x)] = momentY[int(cumX + x)] + float((float(globalForceMatrix[(2 * k)]) * ((x+cumX)-self.SupportLocationArray[k])))

                        if len(pointLoadIndexMoment) > 0:
                            for j in range(0, len(pointLoadIndexMoment)):

                                if x > ((pointLoadLocationArray[pointLoadIndexMoment[j]] * 12)):
                                    momentY[int(cumX + x)] = momentY[int(cumX + x)]+(float(pointLoadMagnitudeArray[pointLoadIndexMoment[j]] * (((x+cumX) - ((pointLoadLocationArray[pointLoadIndexMoment[j]] * 12))) ** 1)))

                        if len(distLoadIndexMoment) > 0:

                            for j in range(0, len(distLoadIndexMoment)):

                                loadLength = float(distributedLoadEndLocationArray[distLoadIndexMoment[j]] - distributedLoadStartLocationArray[distLoadIndexMoment[j]])*12.0

                                if (distributedLoadStartMagnitudeArray[distLoadIndexMoment[j]] == distributedLoadEndMagnitudeArray[distLoadIndexMoment[j]]):

                                    if (x+cumX) > ((((distributedLoadStartLocationArray[distLoadIndexMoment[j]]*12)))):
                                        momentY[int(cumX + x)] = momentY[int(cumX + x)]+(float((distributedLoadStartMagnitudeArray[distLoadIndexMoment[j]]/24)*(((x+cumX)-((distributedLoadStartLocationArray[distLoadIndexMoment[j]]*12)))**2)))
                                    if ((x+cumX) > ((distributedLoadStartLocationArray[distLoadIndexMoment[j]]*12)+loadLength)):
                                        compX = ((distributedLoadStartLocationArray[distLoadIndexMoment[j]]*12)+loadLength)
                                        momentY[int(cumX + x)] = momentY[int(cumX + x)] - (float(( distributedLoadStartMagnitudeArray[distLoadIndexMoment[j]] / 24) * (((x+cumX) - compX) ** 2)))

                                elif (distributedLoadStartMagnitudeArray[distLoadIndexMoment[j]] < distributedLoadEndMagnitudeArray[distLoadIndexMoment[j]]):

                                    if ((distributedLoadStartLocationArray[distLoadIndexMoment[j]]*12)+loadLength-cumX) < spanLength[i]:

                                        compX = ((distributedLoadStartLocationArray[distLoadIndexMoment[j]]*12)+loadLength-cumX)
                                        w = (float((distributedLoadEndMagnitudeArray[distLoadIndexMoment[j]]-distributedLoadStartMagnitudeArray[distLoadIndexMoment[j]]))/12.0)
                                        k = (w/loadLength)

                                        if x > ((distributedLoadStartLocationArray[distLoadIndexMoment[j]] * 12) - cumX):
                                            momentY[int(cumX + x)] = momentY[int(cumX + x)] + (float(((distributedLoadStartMagnitudeArray[distLoadIndexMoment[j]]/24) * ((x-((distributedLoadStartLocationArray[distLoadIndexMoment[j]]*12)-cumX))**2))+((k/6)*((x-((distributedLoadStartLocationArray[distLoadIndexMoment[j]]*12)-cumX))**3))))
                                        if x > ((distributedLoadStartLocationArray[distLoadIndexMoment[j]] * 12) + loadLength - cumX):
                                            momentY[int(cumX + x)] = momentY[int(cumX + x)] - (float(((distributedLoadStartMagnitudeArray[distLoadIndexMoment[j]]/24) * ((x-compX)**2))-(-1*(k/6)*((x-compX)**3)))) - ((w*((x-compX)**2))/2)

                                    else:
                                        if x > ((((distributedLoadStartLocationArray[distLoadIndexMoment[j]] * 12) - cumX))):
                                            k = float((distributedLoadEndMagnitudeArray[distLoadIndexMoment[j]]-distributedLoadStartMagnitudeArray[distLoadIndexMoment[j]])/loadLength)/12.0
                                            momentY[int(cumX + x)] = momentY[int(cumX + x)] + (float(((distributedLoadStartMagnitudeArray[distLoadIndexMoment[j]]/24) * ((x-((distributedLoadStartLocationArray[distLoadIndexMoment[j]]*12)-cumX))**2))+((k/6)*((x-((distributedLoadStartLocationArray[distLoadIndexMoment[j]]*12)-cumX))**3))))

                                elif (distributedLoadStartMagnitudeArray[distLoadIndexMoment[j]] > distributedLoadEndMagnitudeArray[distLoadIndexMoment[j]]):

                                    if ((distributedLoadStartLocationArray[distLoadIndexMoment[j]]*12)+loadLength-cumX) < spanLength[i]:

                                        compX = ((distributedLoadStartLocationArray[distLoadIndexMoment[j]]*12)+loadLength-cumX)
                                        w = (float((distributedLoadStartMagnitudeArray[distLoadIndexMoment[j]]-distributedLoadEndMagnitudeArray[distLoadIndexMoment[j]]))/12.0)
                                        k = (w/loadLength)

                                        if x > ((((distributedLoadStartLocationArray[distLoadIndexMoment[j]] * 12) - cumX))):
                                            momentY[int(cumX + x)] = momentY[int(cumX + x)] + (float(((distributedLoadEndMagnitudeArray[distLoadIndexMoment[j]]/24) * ((x-((distributedLoadStartLocationArray[distLoadIndexMoment[j]]*12)-cumX))**2))+(-1*(k/6)*((x-((distributedLoadStartLocationArray[distLoadIndexMoment[j]]*12)-cumX))**3))+((w/2) * ((x-((distributedLoadStartLocationArray[distLoadIndex[j]]*12)-cumX))**2))))
                                        if x > ((distributedLoadStartLocationArray[distLoadIndexMoment[j]] * 12) + loadLength - cumX):
                                            momentY[int(cumX + x)] = momentY[int(cumX + x)] - (float(((distributedLoadEndMagnitudeArray[distLoadIndexMoment[j]]/24) * ((x-compX)**2))+(-1*(k/6)*(((x)-compX)**3)) + ((w/2)*(((x)-compX)**2))-((w/2) * ((x-compX)**2))))
                                    else:
                                        if x > ((((distributedLoadStartLocationArray[distLoadIndexMoment[j]] * 12) - cumX))):
                                            k = float((distributedLoadStartMagnitudeArray[distLoadIndexMoment[j]]-distributedLoadEndMagnitudeArray[distLoadIndexMoment[j]])/loadLength)/12.0
                                            momentY[int(cumX + x)] = momentY[int(cumX + x)] + (float(((distributedLoadEndMagnitudeArray[distLoadIndexMoment[j]]/24) * ((x-((distributedLoadStartLocationArray[distLoadIndexMoment[j]]*12)-cumX))**2))+((k/6)*((x-((distributedLoadStartLocationArray[distLoadIndexMoment[j]]*12)-cumX))**3))))

                    if x != spanLength[i]:
                        momentX.append(float(cumX + x + 1))
                        momentY.append(0.0)

                if i > 0:
                    cumX = cumX + spanLength[i]
                else:
                    cumX = spanLength[i]

            for i in range(0, len(spanForceMatrix)+1):

                Label(self.master, text="Support " + str(i+1) + " reaction is " + str(float(globalForceMatrix[(2*i)])) + " kips").place(x=350, y=150 + (25*(4*i)), anchor="center")
                Label(self.master, text="Support " + str(i + 1) + " moment is " + str(float(globalForceMatrix[(2 * i)+1])) + " kip*in").place(x=350, y=175 + (25 * (4*i)), anchor="center")

                Label(self.master, text="Support " + str(i+1) + " vertical deflection is " + str(float(dispMatrix[(2*i)])) + " inches").place(x=350, y=200 + (25*(4*i)), anchor="center")
                Label(self.master, text="Support " + str(i+1) + " rotation is " + str(float(dispMatrix[(2*i)+1])) + " radians").place(x=350, y=225 + (25*(4*i)), anchor="center")

            self.ax1.clear()
            self.ax1.plot(shearX, shearY)

            self.ax1.minorticks_on()
            self.ax1.grid(b=True, which='major', color='k', linestyle='-', alpha=0.3)
            self.ax1.grid(b=True, which='minor', color='g', linestyle='-', alpha=0.1)
            self.ax1.set_ylabel('V (kips)')
            self.ax1.set_xlabel('L (in)')

            shearX = np.array(shearX, dtype=float)
            shearY = np.array(shearY, dtype=float)

            self.ax1.fill_between(shearX, np.zeros(len(shearX), dtype=float), shearY, facecolor='red', alpha=0.2)

            self.ax2.clear()
            self.ax2.plot(momentX, momentY)

            self.ax2.minorticks_on()
            self.ax2.grid(b=True, which='major', color='k', linestyle='-', alpha=0.3)
            self.ax2.grid(b=True, which='minor', color='g', linestyle='-', alpha=0.1)
            self.ax2.set_ylabel('M (kips*in)')
            self.ax2.set_xlabel('L (in)')

            maxMoment = max(momentY)
            xMaxMoment = momentY.index(maxMoment)
            xPosMaxMoment = momentX[xMaxMoment]

            minMoment = min(momentY)
            xMinMoment = momentY.index(minMoment)
            xPosMinMoment = momentX[xMinMoment]

            momentX = np.array(momentX, dtype=float)
            momentY = np.array(momentY, dtype=float)
            momentY = momentY.flatten()

            self.ax2.fill_between(momentX, np.zeros(len(momentX), dtype=float), momentY, facecolor='blue', alpha=0.2)

            self.ax2.annotate("Max = " + str(maxMoment), xy=(xPosMaxMoment, maxMoment))
            self.ax2.annotate("Min = " + str(minMoment), xy=(xPosMinMoment, minMoment))

            #print(dispMatrix)
            #print(globalForceMatrix)

        self.after(500, self.Analysis)

main = Tk()
main.geometry("700x500")
app = Application(master=main)
app.Analysis()
app.mainloop()