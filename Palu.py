"""
@author: agus daud
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

# unit conversion
tsf = 105.6
ft = 0.0328084 # cm
m = 0.3048 # ft

class Cohesionless:
    """docstring for Cohesionless."""

    def __init__(self, name="cohesionless"):
        self.name = name

    def fs_ReeseONeil(N,z,po_ef):
        if N >= 15:
            beta = 1.5-(0.245*(z**0.5))
        else:
            beta = (N/15)*1.5-(0.245*(z**0.5))

        if beta <= 0.25:
            beta = 0.25
        elif beta >= 1.2:
            beta = 1.2

        fs = beta*po_ef
        return fs

    def fs_ReeseWright(N):
        if N <= 53:
            fs = (N/34)*tsf
        elif N <= 100:
            fs = (((N-53)/450) + 1.6) * tsf
        else:
            fs = 1.7*tsf
        return fs

    def fs_Meyerhof(N):
        fs = (N/100)*tsf
        return fs

    def fs_QuirosReese(N):
        fs_max = 2*tsf
        fs = (0.026*N)*tsf
        return min(fs,fs_max)

    def fs_ToumaReese(phi,Lp,po_ef):
        phi = np.radians(phi)
        if Lp < 25*ft:
            K = 0.7
        elif Lp < 40*ft:
            K = 0.6
        else:
            K = 0.5

        fs_max = 2.5*tsf
        fs = K*po_ef*np.tan(phi)

        return min(fs,fs_max)

    def fs_API(N,po_ef):
        if N < 10:
            # loose sand
            beta = 0
            fs_max = 0
        elif N < 30:
            # medium sand
            beta = 0.29
            fs_max = 67
        elif N < 50:
            # Dense sand
            beta = 0.37
            fs_max = 81
        elif N >= 50:
            # Very dense sand
            beta = 0.46
            fs_max = 96

        fs =beta*po_ef

        return min(fs,fs_max)

    def fs_Eurocode(phi,po_ef):
        phi = np.radians(phi)
        K = 0.7
        fs = K*po_ef*np.tan(phi)
        return fs

    def fs_FHWA(N60,phi,po_ef):
        phi = np.radians(phi)
        Kp = np.tan(np.radians(45)+(phi/2))**2
        m = 0.8 # assume silty sand to sandy silt
        pa = 100 # kPa
        pc = (0.47*(N60)**m)*100
        # pc = 0.15*N60*100

        beta_max = Kp*np.tan(phi)
        beta = (1 - np.sin(phi))*((pc/po_ef)**np.sin(phi))*np.tan(phi)
        fs_max = po_ef*beta_max
        fs = po_ef*beta

        return min(fs,fs_max)

    def fb_ReeseONeil(N):
        if N > 50:
            N = 50
        fb_max = 2900 # ??? sumber dari mana
        fb = 57.5*N

        return min(fb,fb_max)

    def fb_ReeseWright(N):
        if N <= 60:
            fb = (2/3)*N*tsf
        else:
            fb = 40*tsf
        return fb

    def fb_Meyerhof(N,po_ef,dp):
        N_corr = 0.77*np.log10(20/po_ef)*N # Agus last
        qp = (2*N_corr*dp)/(15*dp)
        qp_max_sand = (4/3)*N_corr # sand
        qp_max_silt = N_corr # silt non plastic

        return min(qp,qp_max_silt)*tsf # tsf

    def fb_ToumaReese(N,dp):
        if dp <= 50:
            K = 0.6
        else:
            K = 1.0

        if N < 10:
            # loose sand
            fb = 0
        elif N <=30 :
            # Medium sand
            fb = (16*tsf)/K
        else:
            # dense
            fb = (40*tsf)/K

        return fb

    def fb_API(N,po_ef):
        if N < 10:
            # loose sand
            Nq = 0
            fb_max = 0
        elif N < 30:
                # medium sand
            Nq = 12
            fb_max = 60
        elif N < 50:
            # Dense sand
            Nq = 20
            fb_max = 100
        elif N >= 50:
            Nq = 40
            fb_max = 200

        fb = Nq*po_ef

        return min(fb,fb_max)

    def fb_Eurocode(phi,po_ef):
        phi = np.radians(phi)
        Nq = round(np.e**(np.pi*np.tan(phi)) * np.tan(np.radians(45) + phi/2)**2,2)
        fb = Nq*po_ef

        return fb

class Cohesive:
    """docstring fo Cohesive."""

    def __init__(self, name="cohesive"):
        self.name = name

    def fb_Skempton(Su):
        Nc = 9
        fb = Nc*Su
        return fb

    def fb_ReeseONeil(Su):
        if Su > 96:
            Nc = 9
        elif Su > 25:
            Nc = 8
        elif Su <= 24:
            Nc = 6.5
        fb = Nc*Su
        return fb

    def fs_Skempton(Su):
        alpha = 0.45
        fs = alpha*Su
        return fs

    def fs_ReeseONeil(Su):
        if Su < 150:
            alpha = 0.55
        else:
            alpha = 0.55-(0.1*((Su/100) - 1.5))
        fs = alpha*Su
        return fs

    def fs_FHWA(Su,po_ef):
        Su_CIUC = Su / (0.911 + (0.499*np.log10(Su/po_ef)))
        pa = 100
        alpha = 0.3 + (0.17/(Su_CIUC/pa))
        fs = alpha*Su_CIUC
        return fs

    def fs_API(Su,po_ef):
        psi = Su/po_ef

        if psi <= 1.0:
            alpha = 0.5*psi**(-0.5)
        else:
            alpha = 0.5*psi**(-0.25)

        alpha_max = 1
        fs_max = Su*alpha_max
        fs = Su*alpha
        return min(fs,fs_max)

class Drilled:
    """docstring for drilled pile."""

    def __init__(self, name="P1", Lp=10, d_cap = 0, D=0.4, GWL=99, gw=9.81, wc=24,ignore_Ltop=0,ignore_Lbot=0,decimal=2,dz=0.5,Rzd=None):
        self.name = name
        self.Lp = Lp
        self.dp = D # cm
        self.GWL = GWL
        self.gw = gw
        self.wc = wc
        self.As = round(3.14*(self.dp/100),3)
        self.Ab = round(0.25*3.14*(self.dp/100)**2,3)
        self.d_cap = d_cap
        self.ignore_Ltop = ignore_Ltop
        self.ignore_Lbot = ignore_Lbot
        self.Q_service = 0
        self.dz = dz
        self.decimal = decimal
        if Rzd != None:
            self.cons_po_z = round(Rzd*(self.dp/100),2)
        else:
            self.cons_po_z = 999

    def insert_row(self,z,data):
        df = data
        index = df.loc[df['z'] >= z].index[0]

        if df.loc[index,'z'] == z:
            df2 = df
        else:
            df.loc[index-0.5] = df.loc[index]
            df.loc[index-0.5,'z'] = z
            df.sort_values(by=['z'])
            df1 = df.sort_index()
            df2 = df1.reset_index(drop=True)
        return df2

    def read_excel(self,data):
        self.raw_df = pd.read_excel(data)
        self.df = pd.read_excel(data)
        self.df = self.insert_row(self.GWL,self.df)
        self.df = self.insert_row(self.Lp,self.df)
        self.df = self.insert_row(self.d_cap,self.df)
        
        if self.ignore_Ltop > 0:
            self.df = self.insert_row(self.ignore_Ltop,self.df)
        if self.ignore_Lbot > 0:
            self.df = self.insert_row(self.Lp-self.ignore_Lbot,self.df)

        list_z = np.arange(0,self.Lp+self.dz,self.dz)
        last_z =  self.df["z"].iloc[-1]
        for z in list_z:
            if z < last_z:
                self.df = self.insert_row(z,self.df)
            else:
                pass

        # constant effectif stress
        if self.cons_po_z < last_z:
            self.df = self.insert_row(self.cons_po_z ,self.df)
        self.df = self.df.reset_index(drop=True)

        # self.insert_row(self.Lp)

    def SPTtoPhi(self,N):
        # phi = ((12*N)**0.5) + 20
        phi = 27.5 + (9.2*np.log10(N))
        return round(phi,2)

    def SPT_undrainedStrength(self,N):
        Su = 6.25*N
        return round(Su,self.decimal)

    def sand_SPT_unitweight(self,N):
        if N < 10:
            g = 18
        else:
            g = 18.5
        return g

    def clay_SPT_unitweight(self,N):
        Su = 6.25*N
        if Su < 25:
            g = 17
        else:
            g = 20
        return g

    def defineParams(self):
        for i, category in enumerate(self.df['category']):
            if category == 'cohesionless':
                phi = self.df.loc[i,'phi']
                if not phi > 0:
                    N = self.df.loc[i,'N']
                    self.df.loc[i,'phi'] = self.SPTtoPhi(N)
                g = self.df.loc[i,'g']
                if not g > 0:
                    N = self.df.loc[i,'N']
                    self.df.loc[i,'g'] = self.sand_SPT_unitweight(N)
            else:
                Su = self.df.loc[i,'Su']
                if not Su > 0:
                    N = self.df.loc[i,'N']
                    self.df.loc[i,'Su'] = self.SPT_undrainedStrength(N)
                g = self.df.loc[i,'g']
                if not g > 0:
                    N = self.df.loc[i,'N']
                    self.df.loc[i,'g'] = self.clay_SPT_unitweight(N)

    def overburden(self):
        self.df.loc[0,'h'] = 0
        self.df.loc[0,'po'] = 0
        self.df.loc[0,'u'] = 0
        self.df.loc[0,'po_ef'] = 0

        for i, g in enumerate(self.df['g']):
            if i > 0:
                po_prev = self.df.loc[i-1,'po']
                z_previous = self.df.loc[i-1,'z']
                z_current = self.df.loc[i,'z']

                h = (z_current - z_previous)
                po = po_prev + (h*g) # total overburden

                if z_current <= self.GWL:
                    u = 0
                elif z_current > self.GWL:
                    u = (z_current - self.GWL)*self.gw # pore water pressure

                self.df.loc[i,'h'] = h
                self.df.loc[i,'po'] = po
                self.df.loc[i,'u'] = u
                self.df.loc[i,'po_ef'] = po - u
            else:
                pass

        if self.cons_po_z < self.df["z"].iloc[-1]:
            i_const_po = self.df.loc[self.df['z'] == self.cons_po_z].index[0]
            self.df.loc[i_const_po:,'po_ef'] = self.df.loc[i_const_po,'po_ef']

    def intrv_params(self):
        i_start = self.df.loc[self.df['z'] == self.d_cap].index[0]
        i_end = self.df.loc[self.df['z'] == self.Lp].index[0]
        mid_df = pd.DataFrame(columns = ['interval' ,'category', 'z', 'h' , 'po_ef', 'Su', 'phi', 'N'])
        for i in range(i_start,i_end):
            # depth and thickness
            z_start = self.df.loc[i,'z']
            z_end = self.df.loc[i+1,'z']
            intrv = f"{round(z_start,2)} - {round(z_end,2)}"
            h = z_end - z_start
            z_mid = z_start + (h/2)
            # overburden
            po_ef_start = self.df.loc[i,'po_ef']
            po_ef_end = self.df.loc[i+1,'po_ef']
            po_ef_avg = round(0.5*(po_ef_end + po_ef_start),2)
            # strength params
            category = self.df.loc[i+1,'category']
            Su = self.df.loc[i+1,'Su']
            phi = self.df.loc[i+1,'phi']
            N = self.df.loc[i+1,'N']

            line = [intrv,category,z_mid,h,po_ef_avg,Su,phi,N]
            mid_df.loc[i] = line

        mid_df = mid_df.reset_index(drop=True)
        self.mid_df = mid_df

    def skin_resistance(self,ESA,TSA):
        null = np.zeros(len(self.mid_df['h']))
        df = self.mid_df
        df['fs'] = null
        df['Qs'] = null

        for i, po_ef in enumerate(self.mid_df['po_ef']):
            category = self.mid_df.loc[i,'category']
            h = self.mid_df.loc[i,'h']

            if category == 'cohesionless':
                z = self.mid_df.loc[i,'z']
                N = self.mid_df.loc[i,'N']
                po_ef = self.mid_df.loc[i,'po_ef']
                phi = self.mid_df.loc[i,'phi']

                if ESA == "Reese and O'Neil":
                    fs = Cohesionless.fs_ReeseONeil(N,z,po_ef)
                elif ESA == "Meyerhof":
                    fs = Cohesionless.fs_Meyerhof(N)
                elif ESA == "Reese and Wright":
                    fs = Cohesionless.fs_ReeseWright(N)
                elif ESA == "Touma and Reese":
                    fs = Cohesionless.fs_ToumaReese(phi,self.Lp,po_ef)
                elif ESA == "Quiros and Reese":
                    fs = Cohesionless.fs_QuirosReese(N)
                elif ESA == "FHWA":
                    fs = Cohesionless.fs_FHWA(N,phi,po_ef)
                elif ESA == "Eurocode":
                    fs = Cohesionless.fs_Eurocode(phi,po_ef)
                elif ESA == "API":
                    fs = Cohesionless.fs_API(N,po_ef)
                # Skin resistance
                Qs = fs*self.As*h
            elif category == 'cohesive':
                Su = self.mid_df.loc[i,'Su']
                po_ef = self.mid_df.loc[i,'po_ef']

                if TSA == "Skempton":
                    fs = Cohesive.fs_Skempton(Su)
                elif TSA == "Reese and O'Neil":
                    fs = Cohesive.fs_ReeseONeil(Su)
                elif TSA == "FHWA":
                    fs = Cohesive.fs_FHWA(Su,po_ef)
                elif TSA == "API":
                    fs = Cohesive.fs_API(Su,po_ef)
                Qs = fs*self.As*h

            df.loc[i,'fs'] = round(fs,self.decimal)
            df.loc[i,'Qs'] = round(Qs,self.decimal)

        if self.ignore_Ltop > 0:
            index = self.df.loc[self.df['z'] == self.ignore_Ltop].index[0] - 1
            df.loc[:index,"fs"] = 0
            df.loc[:index,"Qs"] = 0

        if self.ignore_Lbot > 0:
            index = self.df.loc[self.df['z'] == self.Lp-self.ignore_Lbot].index[0]
            df.loc[index:,"fs"] = 0
            df.loc[index:,"Qs"] = 0

        Qs_tot = round(sum(df['Qs']),self.decimal)
        return Qs_tot, df

    def end_bearing(self,ESA,TSA):
        index = self.df.loc[self.df['z'] == self.Lp].index[0]
        end_soil = self.df.loc[index,'category']

        if end_soil == 'cohesionless':
            N = self.df.loc[index,'N']
            po_ef = self.df.loc[index,'po_ef']/tsf
            dp = self.dp
            if ESA == "Meyerhof":
                fb = Cohesionless.fb_Meyerhof(N,po_ef,dp)
            elif ESA == "Reese and O'Neil":
                fb = Cohesionless.fb_ReeseONeil(N)
            elif ESA == "Reese and Wright":
                fb = Cohesionless.fb_ReeseWright(N)
            elif ESA == "Touma and Reese":
                fb = Cohesionless.fb_ToumaReese(N,dp)
            elif ESA == "API":
                fb = Cohesionless.fb_API(N,po_ef)
            elif ESA == "Eurocode":
                fb = Cohesionless.fb_Eurocode(phi,po_ef)
            Qb = fb*self.Ab
        elif end_soil == 'cohesive':
            Su = self.df.loc[index,'Su']
            if TSA == "Skempton":
                fb = Cohesive.fb_Skempton(Su)
            elif TSA == "Reese and O'Neil":
                fb = Cohesive.fb_ReeseONeil(Su)
            Qb = fb*self.Ab

        return round(Qb,self.decimal)

    def limitStateDisplacement(self,S=2.5):
        delta_B = round((S/self.dp)*100,self.decimal) # normalized torable displacement
        x = delta_B
        category = self.mid_df.iloc[-1]['category']
        if category == "cohesive":
            if delta_B <= 0.4:
                norm_axial = x * (50/0.4)
            elif delta_B <= 4:
                norm_axial = (0.2688*x**5) - (3.7324*x**4) + (20.135*x**3) - (55.318*x**2) + (88.621*x) + 22.205
            elif delta_B > 4:
                norm_axial = 100
        elif category == "cohesionless":
            if delta_B <= 0.4:
                norm_axial = x*(50/0.4)
            elif delta_B <= 4:
                norm_axial = (0.2717*x**5) - (3.6614*x**4) + (19.022*x**3) - (47.779*x**2) + (69.288*x) + 28.949
            else:
                norm_axial = 100 + (x-4)*(59/6)
        return round(norm_axial,self.decimal), delta_B

    def services_load(self,S=2.5):
        self.S_serv = S
        self.norm_axial, self.delta_B = self.limitStateDisplacement(S=S)
        self.Q_threshold = self.Qult
        self.Q_service = round((self.Q_threshold*(self.norm_axial/100)) - self.Wp,0)
        if self.Q_service < 0:
            self.Q_service = 0

    def TZ_curve(self,S=2.5,test_data=None):
        list_settlement = np.arange(0,S,0.01)
        list_Q = []
        for S in list_settlement:
            self.services_load(S)
            list_Q.append(self.Q_service)
        df = pd.DataFrame(dict(Q = list_Q, S = list_settlement))
        # f, ax = plt.subplots()
        # ax.plot(list_Q,list_settlement)
        # ax.set_xlabel("Load Q (kN)")
        # ax.set_ylabel("Settlement (cm)")
        # ax.invert_yaxis()

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=list_Q,y=list_settlement,name = 'Estimate',connectgaps=True))
        if test_data != None:
            x = test_data[0]
            y = test_data[1]
            fig.add_trace(go.Scatter(x=x,y=y,name = 'Pile Test',connectgaps=True))

        # fig = px.line(df,x='Q', y='S', title=f'TZ curve',width=520,height=600)
        # fig = px.line(df2,x='Q', y='S', title=f'TZ curve',width=520,height=600)
        fig.update_layout(
            width=600,
            height=600,
            title="Load - Settlement Curve",
            yaxis_title="Settlement (cm)",
            xaxis_title="Load (kN)")
        fig['layout']['yaxis']['autorange'] = "reversed"
        fig.show()

    def expected_settl(self,Q):
        self.Qu = Q + self.Wp
        self.Qserv_expected = Q
        self.norm_axial_expected = round((self.Qu/self.Q_threshold)*100,self.decimal)

        S_init = 0.01
        while True:
            norm_axial, delta_B = self.limitStateDisplacement(S=S_init)
            if abs(self.norm_axial_expected-norm_axial) < 0.5:
                self.Sserv_expected = round(S_init,2)
                self.delta_B_expected = delta_B
                break
            else:
                S_init += 0.001

    def solve(self,FS=2.5,Qb = ["Reese and Wright","Skempton"],Qs=["Meyerhof","Reese and O'Neil"],ignore_fs_end=False):
        self.ignore_fs_end = ignore_fs_end
        self.FS = FS

        self.defineParams()
        self.overburden()
        self.intrv_params()

        # water_pressure = self.df.iloc[-1]['u']
        self.Wp = round(self.Ab*self.Lp*self.wc) - round(self.Ab*(self.Lp-self.GWL)*9.81)

        method1, method2 = Qs
        self.Qs, self.table = self.skin_resistance(ESA=method1, TSA=method2)

        method1, method2 = Qb
        self.Qb = self.end_bearing(ESA=method1, TSA=method2)

        self.Qup = (0.75 * round(self.Qs + self.Wp,2))/self.FS

        self.Qult = round(self.Qs + self.Qb,2)
        self.Qa = round(self.Qult/self.FS,2)

    def Msolve(self,FS=2.5,methods=None,ignore_fs_end=False,S=0):
        self.FS = FS
        self.defineParams()
        self.overburden()
        self.intrv_params()
        if methods != None:
            self.methods = methods
        else:
            # Skin friction
            ESA_fr = ["Meyerhof","Reese and Wright","Reese and O'Neil","Touma and Reese","Quiros and Reese","FHWA","API"]
            TSA_fr = ["Reese and O'Neil","Reese and O'Neil","Reese and O'Neil","Reese and O'Neil","Reese and O'Neil","FHWA","API"]
            # End bearing
            ESA_end = ["Meyerhof","Reese and Wright","Reese and O'Neil","Touma and Reese","Reese and O'Neil","Reese and O'Neil","API"]
            TSA_end = ["Reese and O'Neil","Reese and O'Neil","Reese and O'Neil","Reese and O'Neil","Reese and O'Neil","Reese and O'Neil","Skempton"]

            self.methods = ESA_fr

        self.list_Qs = []
        self.list_Qb = []
        self.list_Qult = []
        self.list_Qa = []
        self.list_table = []
        self.list_Q_service = []
        self.list_Qup = []

        for i in range(7):
            self.solve(FS=FS,Qb=[ESA_end[i],TSA_end[i]],Qs=[ESA_fr[i],TSA_fr[i]],ignore_fs_end=ignore_fs_end)
            if S>0:
                self.services_load(S=S)
                self.list_Q_service.append(self.Q_service)
            self.list_Qs.append(self.Qs)
            self.list_Qb.append(self.Qb)
            self.list_Qult.append(self.Qult)
            self.list_Qa.append(self.Qa)
            self.list_table.append(self.table)
            self.list_Qup.append(self.Qup)

    def plot_result(self,param="Qu"):
        if param == "Qult":
            y = self.list_Qult
        elif param == "Qa":
            y = self.list_Qa
        elif param == "Qs":
            y = self.list_Qs
        elif param == "Qb":
            y = self.list_Qb
        elif param == "Q_service":
            y = self.list_Q_service
        elif param == "Q_uplift":
            y = self.list_Qup
        fig = px.bar(x=self.methods, y=y, text_auto='.3s',title="Pile Capacity", labels={'x': "Methods", 'y':f"{param} (kN)"},width=600, height=500)
        fig.update_traces(width=0.4,textfont_size=12, textangle=0, textposition="outside", cliponaxis=False)
        fig.show()

    def info(self):
        print(f'Qs = {self.Qs} kN')
        print(f'Qb = {self.Qb} kN')
        print(f'Wp = {round(self.Wp,2)} kN')
        print(f'Qult = {self.Qult} kN')
        print(f'Qa = {self.Qa} kN')
        print(f'Q uplift = {round(self.Qup,0)} kN')

    def service_condition(self):
        print(f"Norm. axial force = {self.norm_axial} %")
        print(f"Norm. displacement = {self.delta_B} %")
        print(f"Failure threshold = {self.Q_threshold} kN")
        print(f"Services load, Q = {self.Q_service} kN")
        print(f"Settlement = {self.S_serv} cm")

    def expected_condition(self):
        print(f"Norm. axial force = {self.norm_axial_expected} %")
        print(f"Norm. displacement = {self.delta_B_expected} %")
        print(f"Failure threshold = {self.Qult} kN")
        print(f"Services load, Q = {self.Qserv_expected} kN")
        print(f"Settlement expected = {self.Sserv_expected} cm")

    def plot(self,param="fs"):
        fig = px.line(self.table, x=param, y="z", title=f'Graph: depth z (m) vs {param} kPa',width=320,height=600)
        # fig.update_traces(width=2,height=6)
        fig['layout']['yaxis']['autorange'] = "reversed"
        fig.show()

    def to_csv(self,name):
        self.table.to_csv(name)

class Driven_Cohesionless:
    """docstring for Cohesionless."""

    def __init__(self, name="cohesionless"):
        self.name = name

    def fb_Meyerhof(N60,Lp,dp):
        pa = 100
        fb = 0.4*N60*(Lp/dp)*pa
        fb_max = 4*N60*pa
        return min(fb,fbmax)

    def fs_Meyerhof(N60):
        pa = 100
        fs = (1/50)*pa*N60
        return fs

class Driven_Cohesive:
    """docstring for Cohesionless."""

    def __init__(self, name="cohesionless"):
        self.name = name

    def fb_Skempton(Su):
        Nc = 9
        fb = 9*Su
        return fb

    def fs_API(Su):
        if Su <= 25:
            alpha = 1
        elif Su >= 70:
            alpha = 0.5
        else:
            alpha = 1 - (Su-25)/90
        fs = alpha*Su
        return fs

class Driven:
    """docstring for drilled pile."""

    def __init__(self, name="P1", Lp=10, d_cap = 0, D=0.4, GWL=99, gw=9.81, wc=24,ignore_Ltop=0,ignore_Lbot=0,decimal=2,dz=0.5,Rzd=None):
        self.name = name
        self.Lp = Lp
        self.dp = D # cm
        self.GWL = GWL
        self.gw = gw
        self.wc = wc
        self.As = round(3.14*(self.dp/100),3)
        self.Ab = round(0.25*3.14*(self.dp/100)**2,3)
        self.d_cap = d_cap
        self.ignore_Ltop = ignore_Ltop
        self.ignore_Lbot = ignore_Lbot
        self.Q_service = 0
        self.dz = dz
        self.decimal = decimal
        if Rzd != None:
            self.cons_po_z = round(Rzd*(self.dp/100),2)
        else:
            self.cons_po_z = 999

    def insert_row(self,z,data):
        df = data
        index = df.loc[df['z'] >= z].index[0]

        if df.loc[index,'z'] == z:
            df2 = df
        else:
            df.loc[index-0.5] = df.loc[index]
            df.loc[index-0.5,'z'] = z
            df.sort_values(by=['z'])
            df1 = df.sort_index()
            df2 = df1.reset_index(drop=True)
        return df2

    def read_excel(self,data):
        self.raw_df = pd.read_excel(data)
        self.df = pd.read_excel(data)
        self.df = self.insert_row(self.GWL,self.df)
        self.df = self.insert_row(self.Lp,self.df)
        self.df = self.insert_row(self.d_cap,self.df)
        if self.ignore_Ltop > 0:
            self.df = self.insert_row(self.ignore_Ltop,self.df)
        if self.ignore_Lbot > 0:
            self.df = self.insert_row(self.Lp-self.ignore_Lbot,self.df)

        list_z = np.arange(0,self.Lp+self.dz,self.dz)
        last_z =  self.df["z"].iloc[-1]
        for z in list_z:
            if z < last_z:
                self.df = self.insert_row(z,self.df)
            else:
                pass

        # constant effectif stress
        if self.cons_po_z < last_z:
            self.df = self.insert_row(self.cons_po_z ,self.df)
        self.df = self.df.reset_index(drop=True)

        # self.insert_row(self.Lp)

    def SPTtoPhi(self,N):
        # phi = ((12*N)**0.5) + 20
        phi = 27.5 + (9.2*np.log10(N))
        return round(phi,2)

    def SPT_undrainedStrength(self,N):
        Su = 6.25*N
        return round(Su,self.decimal)

    def sand_SPT_unitweight(self,N):
        if N < 10:
            g = 18
        else:
            g = 18.5
        return g

    def clay_SPT_unitweight(self,N):
        Su = 6.25*N
        if Su < 25:
            g = 17
        else:
            g = 20
        return g

    def defineParams(self):
        for i, category in enumerate(self.df['category']):
            if category == 'cohesionless':
                phi = self.df.loc[i,'phi']
                if not phi > 0:
                    N = self.df.loc[i,'N']
                    self.df.loc[i,'phi'] = self.SPTtoPhi(N)
                g = self.df.loc[i,'g']
                if not g > 0:
                    N = self.df.loc[i,'N']
                    self.df.loc[i,'g'] = self.sand_SPT_unitweight(N)
            else:
                Su = self.df.loc[i,'Su']
                if not Su > 0:
                    N = self.df.loc[i,'N']
                    self.df.loc[i,'Su'] = self.SPT_undrainedStrength(N)
                g = self.df.loc[i,'g']
                if not g > 0:
                    N = self.df.loc[i,'N']
                    self.df.loc[i,'g'] = self.clay_SPT_unitweight(N)

    def overburden(self):
        self.df.loc[0,'h'] = 0
        self.df.loc[0,'po'] = 0
        self.df.loc[0,'u'] = 0
        self.df.loc[0,'po_ef'] = 0

        for i, g in enumerate(self.df['g']):
            if i > 0:
                po_prev = self.df.loc[i-1,'po']
                z_previous = self.df.loc[i-1,'z']
                z_current = self.df.loc[i,'z']

                h = (z_current - z_previous)
                po = po_prev + (h*g) # total overburden

                if z_current <= self.GWL:
                    u = 0
                elif z_current > self.GWL:
                    u = (z_current - self.GWL)*self.gw # pore water pressure

                self.df.loc[i,'h'] = h
                self.df.loc[i,'po'] = po
                self.df.loc[i,'u'] = u
                self.df.loc[i,'po_ef'] = po - u
            else:
                pass

        if self.cons_po_z < self.df["z"].iloc[-1]:
            i_const_po = self.df.loc[self.df['z'] == self.cons_po_z].index[0]
            self.df.loc[i_const_po:,'po_ef'] = self.df.loc[i_const_po,'po_ef']

    def intrv_params(self):
        i_start = self.df.loc[self.df['z'] == self.d_cap].index[0]
        i_end = self.df.loc[self.df['z'] == self.Lp].index[0]
        mid_df = pd.DataFrame(columns = ['interval' ,'category', 'z', 'h' , 'po_ef', 'Su', 'phi', 'N'])
        for i in range(i_start,i_end):
            # depth and thickness
            z_start = self.df.loc[i,'z']
            z_end = self.df.loc[i+1,'z']
            intrv = f"{round(z_start,2)} - {round(z_end,2)}"
            h = z_end - z_start
            z_mid = z_start + (h/2)
            # overburden
            po_ef_start = self.df.loc[i,'po_ef']
            po_ef_end = self.df.loc[i+1,'po_ef']
            po_ef_avg = round(0.5*(po_ef_end + po_ef_start),2)
            # strength params
            category = self.df.loc[i+1,'category']
            Su = self.df.loc[i+1,'Su']
            phi = self.df.loc[i+1,'phi']
            N = self.df.loc[i+1,'N']

            line = [intrv,category,z_mid,h,po_ef_avg,Su,phi,N]
            mid_df.loc[i] = line

        mid_df = mid_df.reset_index(drop=True)
        self.mid_df = mid_df

    def skin_resistance(self,ESA,TSA):
        null = np.zeros(len(self.mid_df['h']))
        df = self.mid_df
        df['fs'] = null
        df['Qs'] = null

        for i, po_ef in enumerate(self.mid_df['po_ef']):
            category = self.mid_df.loc[i,'category']
            h = self.mid_df.loc[i,'h']

            if category == 'cohesionless':
                z = self.mid_df.loc[i,'z']
                N = self.mid_df.loc[i,'N']
                po_ef = self.mid_df.loc[i,'po_ef']
                phi = self.mid_df.loc[i,'phi']

                if ESA == "Meyerhof":
                    fs = Driven_Cohesionless.fs_Meyerhof(N)
                # Skin resistance
                Qs = fs*self.As*h
            elif category == 'cohesive':
                Su = self.mid_df.loc[i,'Su']
                po_ef = self.mid_df.loc[i,'po_ef']

                if TSA == "API":
                    fs = Driven_Cohesive.fs_API(Su)
                Qs = fs*self.As*h

            df.loc[i,'fs'] = round(fs,self.decimal)
            df.loc[i,'Qs'] = round(Qs,self.decimal)

        if self.ignore_Ltop > 0:
            index = self.df.loc[self.df['z'] == self.ignore_Ltop].index[0] - 1
            df.loc[:index,"fs"] = 0
            df.loc[:index,"Qs"] = 0

        if self.ignore_Lbot > 0:
            index = self.df.loc[self.df['z'] == self.Lp-self.ignore_Lbot].index[0]
            df.loc[index:,"fs"] = 0
            df.loc[index:,"Qs"] = 0

        Qs_tot = round(sum(df['Qs']),self.decimal)
        return Qs_tot, df

    def end_bearing(self,ESA,TSA):
        index = self.df.loc[self.df['z'] == self.Lp].index[0]
        end_soil = self.df.loc[index,'category']

        if end_soil == 'cohesionless':
            N = self.df.loc[index,'N']
            dp = self.dp
            Lp = self.Lp
            if ESA == "Meyerhof":
                fb = Cohesionless.fb_Meyerhof(N,Lp,dp)
            Qb = fb*self.Ab
        elif end_soil == 'cohesive':
            Su = self.df.loc[index,'Su']
            if TSA == "Skempton":
                fb = Driven_Cohesive.fb_Skempton(Su)
            Qb = fb*self.Ab
        return round(Qb,self.decimal)

    def limitStateDisplacement(self,S=2.5):
        delta_B = round((S/self.dp)*100,self.decimal) # normalized torable displacement
        x = delta_B
        category = self.mid_df.iloc[-1]['category']
        if category == "cohesive":
            if delta_B <= 0.4:
                norm_axial = x * (50/0.4)
            elif delta_B <= 4:
                norm_axial = (0.2688*x**5) - (3.7324*x**4) + (20.135*x**3) - (55.318*x**2) + (88.621*x) + 22.205
            elif delta_B > 4:
                norm_axial = 100
        elif category == "cohesionless":
            if delta_B <= 0.4:
                norm_axial = x*(50/0.4)
            elif delta_B <= 4:
                norm_axial = (0.2717*x**5) - (3.6614*x**4) + (19.022*x**3) - (47.779*x**2) + (69.288*x) + 28.949
            else:
                norm_axial = 100 + (x-4)*(59/6)
        return round(norm_axial,self.decimal), delta_B

    def services_load(self,S=2.5):
        self.S_serv = S
        self.norm_axial, self.delta_B = self.limitStateDisplacement(S=S)
        self.Q_threshold = self.Qult
        self.Q_service = round((self.Q_threshold*(self.norm_axial/100)) - self.Wp,0)
        if self.Q_service < 0:
            self.Q_service = 0

    def TZ_curve(self,S=2.5,test_data=None):
        list_settlement = np.arange(0,S,0.01)
        list_Q = []
        for S in list_settlement:
            self.services_load(S)
            list_Q.append(self.Q_service)
        df = pd.DataFrame(dict(Q = list_Q, S = list_settlement))
        # f, ax = plt.subplots()
        # ax.plot(list_Q,list_settlement)
        # ax.set_xlabel("Load Q (kN)")
        # ax.set_ylabel("Settlement (cm)")
        # ax.invert_yaxis()

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=list_Q,y=list_settlement,name = 'Estimate',connectgaps=True))
        if test_data != None:
            x = test_data[0]
            y = test_data[1]
            fig.add_trace(go.Scatter(x=x,y=y,name = 'Pile Test',connectgaps=True))

        # fig = px.line(df,x='Q', y='S', title=f'TZ curve',width=520,height=600)
        # fig = px.line(df2,x='Q', y='S', title=f'TZ curve',width=520,height=600)
        fig.update_layout(
            width=600,
            height=600,
            title="Load - Settlement Curve",
            yaxis_title="Settlement (cm)",
            xaxis_title="Load (kN)")
        fig['layout']['yaxis']['autorange'] = "reversed"
        fig.show()

    def expected_settl(self,Q):
        self.Qu = Q + self.Wp
        self.Qserv_expected = Q
        self.norm_axial_expected = round((self.Qu/self.Q_threshold)*100,self.decimal)

        S_init = 0.01
        while True:
            norm_axial, delta_B = self.limitStateDisplacement(S=S_init)
            if abs(self.norm_axial_expected-norm_axial) < 0.5:
                self.Sserv_expected = round(S_init,2)
                self.delta_B_expected = delta_B
                break
            else:
                S_init += 0.001

    def solve(self,FS=2.5,Qb = ["Meyerhof","Skempton"],Qs=["Meyerhof","API"],ignore_fs_end=False):
        self.ignore_fs_end = ignore_fs_end
        self.FS = FS

        self.defineParams()
        self.overburden()
        self.intrv_params()

        # water_pressure = self.df.iloc[-1]['u']
        self.Wp = round(self.Ab*self.Lp*self.wc) - round(self.Ab*(self.Lp-self.GWL)*9.81)

        method1, method2 = Qs
        self.Qs, self.table = self.skin_resistance(ESA=method1, TSA=method2)

        method1, method2 = Qb
        self.Qb = self.end_bearing(ESA=method1, TSA=method2)

        self.Qup = (0.75 * round(self.Qs + self.Wp,2))/self.FS

        self.Qult = round(self.Qs + self.Qb,2)
        self.Qa = round(self.Qult/self.FS,2)

    def Msolve(self,FS=2.5,methods=None,ignore_fs_end=False,S=0):
        self.FS = FS
        self.defineParams()
        self.overburden()
        self.intrv_params()
        if methods != None:
            self.methods = methods
        else:
            # Skin friction
            ESA_fr = ["Meyerhof","Reese and Wright","Reese and O'Neil","Touma and Reese","Quiros and Reese","FHWA","API"]
            TSA_fr = ["Reese and O'Neil","Reese and O'Neil","Reese and O'Neil","Reese and O'Neil","Reese and O'Neil","FHWA","API"]
            # End bearing
            ESA_end = ["Meyerhof","Reese and Wright","Reese and O'Neil","Touma and Reese","Reese and O'Neil","Reese and O'Neil","API"]
            TSA_end = ["Reese and O'Neil","Reese and O'Neil","Reese and O'Neil","Reese and O'Neil","Reese and O'Neil","Reese and O'Neil","Skempton"]

            self.methods = ESA_fr

        self.list_Qs = []
        self.list_Qb = []
        self.list_Qult = []
        self.list_Qa = []
        self.list_table = []
        self.list_Q_service = []
        self.list_Qup = []

        for i in range(7):
            self.solve(FS=FS,Qb=[ESA_end[i],TSA_end[i]],Qs=[ESA_fr[i],TSA_fr[i]],ignore_fs_end=ignore_fs_end)
            if S>0:
                self.services_load(S=S)
                self.list_Q_service.append(self.Q_service)
            self.list_Qs.append(self.Qs)
            self.list_Qb.append(self.Qb)
            self.list_Qult.append(self.Qult)
            self.list_Qa.append(self.Qa)
            self.list_table.append(self.table)
            self.list_Qup.append(self.Qup)

    def plot_result(self,param="Qu"):
        if param == "Qult":
            y = self.list_Qult
        elif param == "Qa":
            y = self.list_Qa
        elif param == "Qs":
            y = self.list_Qs
        elif param == "Qb":
            y = self.list_Qb
        elif param == "Q_service":
            y = self.list_Q_service
        elif param == "Q_uplift":
            y = self.list_Qup
        fig = px.bar(x=self.methods, y=y, text_auto='.3s',title="Pile Capacity", labels={'x': "Methods", 'y':f"{param} (kN)"},width=600, height=500)
        fig.update_traces(width=0.4,textfont_size=12, textangle=0, textposition="outside", cliponaxis=False)
        fig.show()

    def info(self):
        print(f'Qs = {self.Qs} kN')
        print(f'Qb = {self.Qb} kN')
        print(f'Wp = {round(self.Wp,2)} kN')
        print(f'Qult = {self.Qult} kN')
        print(f'Qa = {self.Qa} kN')
        print(f'Q uplift = {round(self.Qup,0)} kN')

    def service_condition(self):
        print(f"Norm. axial force = {self.norm_axial} %")
        print(f"Norm. displacement = {self.delta_B} %")
        print(f"Failure threshold = {self.Q_threshold} kN")
        print(f"Services load, Q = {self.Q_service} kN")
        print(f"Settlement = {self.S_serv} cm")

    def expected_condition(self):
        print(f"Norm. axial force = {self.norm_axial_expected} %")
        print(f"Norm. displacement = {self.delta_B_expected} %")
        print(f"Failure threshold = {self.Qult} kN")
        print(f"Services load, Q = {self.Qserv_expected} kN")
        print(f"Settlement expected = {self.Sserv_expected} cm")

    def plot(self,param="fs"):
        fig = px.line(self.table, x=param, y="z", title=f'Graph: depth z (m) vs {param} kPa',width=320,height=600)
        # fig.update_traces(width=2,height=6)
        fig['layout']['yaxis']['autorange'] = "reversed"
        fig.show()

    def to_csv(self,name):
        self.table.to_csv(name)
