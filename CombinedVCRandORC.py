- 👋 Hi, I’m @vishalkoshti
- 👀 I’m interested in ...
- 🌱 I’m currently learning ...
- 💞️ I’m looking to collaborate on ...
- 📫 How to reach me ...

<!---
vishalkoshti/vishalkoshti is a ✨ special ✨ repository because its `README.md` (this file) appears on your GitHub profile.
You can click the Preview link to take a look at your changes.
--->
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from tabulate import tabulate
import CoolProp
from CoolProp.CoolProp import PropsSI
import numpy
import decimal
import math

def work(fluid,mdot,T1,T4,T2,T6):
    T1=round(T1+273.15,4)
    T4=round(T4+273.15,4)
    T2=round(T2+273.15,4)
    T6=round(T6+273.15,4)
    T7=T1
    T5=T4
    #pressures
    P1=round(PropsSI('P','T',T1,'Q',1,fluid)/1e5,4)
    P4=round(PropsSI('P','T',T4,'Q',1,fluid)/1e5,4)
    P6=P5=P3=P4
    P7=P2=P1
    Cg1=round(PropsSI('C','T',T1,'Q',1,fluid),4)
    Cg4=round(PropsSI('C','T',T4,'Q',1,fluid),4)
    Cl5=round(PropsSI('C','T',T5,'Q',0,fluid),4)
    S1=round(PropsSI('S','T',T1,'Q',1,fluid),4)
    S4=round(PropsSI('S','T',T4,'Q',1,fluid),4)
    T3=(S1-S4+(Cg1*(math.log(T2/T1))))/Cg4
    T3=round((math.e**T3)*T4,4)
    S2=round(S1+(Cg1*(math.log(T2/T1,math.e))),4)
    S3=round(S4+(Cg1*(math.log(T3/T4,math.e))),4)
    S5=round(PropsSI('S','T',T5,'Q',0,fluid),4)
    S6=round(PropsSI('S','T',T6,'Q',0,fluid),4)
    H1=round(PropsSI('H','T',T1,'Q',1,fluid),4)
    H4=round(PropsSI('H','T',T4,'Q',1,fluid),4)
    H5=round(PropsSI('H','T',T5,'Q',0,fluid),4)
    H2=round(H1+(Cg1*(T2-T1)),4)
    H3=round(H4+(Cg4*(T3-T4)),4)
    H6=round(H5-(Cl5*(T5-T6)),4)
    H7=H6
    #state4
    H7f=round(PropsSI('H','T',T7,'Q',0,fluid),4)
    H7fg=round(PropsSI('H','T',T7,'Q',1,fluid)-H7f,4)
    x=round((H7-H7f)/H7fg,2)
    S7f=round(PropsSI('S','T',T7,'Q',0,fluid),4)
    S7fg=round(PropsSI('S','T',T7,'Q',1,fluid)-S7f,4)
    S7=round(S7f+(x*S7fg),4)
    COPideal=round(T1/(T4-T1),4)
    COPactual=round((H2-H7)/(H3-H2),4)
    RC=round(((H1-H7)*(mdot)/3500),4)#H1
    CondesorHeat=round(mdot*(H3-H6),4)#check
    Compressorpower=round(mdot*(H3-H2),4)
    #Qin=CondesorHeat+((Cl5*(T5-T6))*mdot)
    #Qin=CondesorHeat+((H5-H6)*mdot)
    Qin=CondesorHeat
    #print(Qin)
    print('**********************Organic Rankine Cycle********************')
    fluid2=input('Enter the fluid for ORC\n')
    print("#########")
    print("NOTE:",PropsSI("Pcrit",fluid2)/1e5,"is the maximum Pressure Allowed")
    print("#########")
    #Qin=float(input('Enter the Q rejected from VCR in KJ/Kg\n'))
    Pa=float(input('Enter the Pressure of Boiler in Bar\n'))
    Pf=float(input('Enter the Condensor Pressure in Bar(Exit of Turbine)\n'))
    Np=float(input('Enter the Efficiency of Pump in %\n'))
    Nt=float(input('Enter the Efficiency of Turbine in %\n'))
    Np=Np/100
    Nt=Nt/100         
    Pa=Pa*100000
    Pf=Pf*100000
    Pb=Pc=Pd=Pa
    Pe=Pg=Pf
    Tsat=round(PropsSI('T','P',Pa,'Q',1,fluid2),4)
    Tf=round(PropsSI('T','P',Pf,'Q',1,fluid2),4)
    Tg=Tf
    Tb=Tc=Tsat
    Hb=round(PropsSI('H','P',Pc,'Q',0,fluid2),4)
    Sb=round(PropsSI('S','P',Pc,'Q',0,fluid2),4)
    Hc=round(PropsSI('H','P',Pc,'Q',1,fluid2),4)
    Sc=round(PropsSI('S','P',Pc,'Q',1,fluid2),4)
    Vg=round(1/(PropsSI('D','P',Pg,'Q',0,fluid2)),4)
    Sf0=round(PropsSI('S','P',Pf,'Q',0,fluid2),4)
    Sf1=round(PropsSI('S','P',Pg,'Q',1,fluid2),4)
    Hg=round(PropsSI('H','P',Pg,'Q',0,fluid2),4)
    Sg=round(PropsSI('S','P',Pg,'Q',0,fluid2),4)
    mdot_ORC=int(1)
    Wp=(Vg*(Pa-Pf))*mdot_ORC
    Ha=(Wp/mdot_ORC)+Hg
    Hd=(Qin/mdot_ORC)+Ha
    #print('ha',Ha)
    #print(Hd)
    #print(Qin/mdot)
    def ORC_sat(Tc,Td,Hc,Hd,Sc,Sd,Se,Sf0,Sf1,mdot_ORC):
        x=round((Se-Sf0)/(Sf1-Sf0),2)
        print('X is ',x)
        def condensorORC(fluid,fluid2,Pa,Pb,Pc,Pd,Pe,Pf,Pg,Ha,Hb,Hc,Hd,He,Hf,Hg,Tb,Tc,Td,Te,Tf,Tg,Sb,Sc,Sd,Se,Sf,T1,T2,T3,T4,T5,T6,T7,S1,S2,S3,S4,S5,S6,S7,P1,P2,P3,P4,P5,P6,P7,H1,H2,H3,H4,H5,H6,H7,COPideal,COPactual,RC,CondesorHeat,Compressorpower,mdot):
            Wt=(Hd-He)*mdot_ORC
            Wp=(Vg*(Pa-Pf))*mdot_ORC
            ORC_Efficiency=round((((Wt-Wp)/(Qin))*100),2)
            SteamRate=round((3600*1000*mdot_ORC)/(Wt-Wp),4)
            Wpactual=Wp/Np
            Wtactual=Wt*Nt
            ORC_Efficiency_actual=round((((Wtactual-Wpactual)/(Qin))*100),2)
            H_a=(Wpactual/mdot_ORC)+Hg
            H_e=Hd-(Wtactual/mdot_ORC)
            Sa=Sg#check
            S_a=round(PropsSI('S','P',Pa,'H',H_a,fluid2),4)
            S_e=round(PropsSI('S','P',Pe,'H',H_e,fluid2),4)
            Ta=round(PropsSI('T','P',Pa,'H',Ha,fluid2),4)#round(PropsSI('T','P',Pa,'H',Ha,fluid2),4)
            #print('Sa',Se)
            #print('S_a',S_e)
    #checking heat exchanger design
            if T4<Ta:
                print('The temperature in Condensor is',T4)
                print('The temperature in Condensor is',Ta)
                print('Heat Exchanger Design not Possible,Error!!!')
            else:
                T_diff_heatexchanger=T4-Ta
                print("****************************************************")
                print('Heat Exchanger Desing possible')
                print('The Temperature Difference is',T_diff_heatexchanger)
                print("****************************************************")
                print('--------------------------------Vapor Compression Cycle--------------------------------')
                table1=[['1',T1,round(S1/1000,4),round(H1/1000,4),P1],
                        ['2',T2,round(S2/1000,4),round(H2/1000,4),P2],
                        ['3',T3,round(S3/1000,4),round(H3/1000,4),P3],
                        ['4',T4,round(S4/1000,4),round(H4/1000,4),P4],
                        ['5',T5,round(S5/1000,4),round(H5/1000,4),P5],
                        ['6',T6,round(S6/1000,4),round(H6/1000,4),P6],
                        ['7',T7,round(S7/1000,4),round(H7/1000,4),P7]]
                print(tabulate(table1,headers=['State Point','Temperature(in K)','Entropy(S)(in KJ/Kg)','Enthaply(H)(KJ/Kg)','Pressure(Bar)']))
                table2=[['Carnot COP',COPideal,'--'],
                        ['Actual COP',COPactual,'--'],
                        ['Refriferation Capacity',RC,'tons'],
                        ['Heat Rejected in Condensor',round(CondesorHeat/1000,4),'KW'],
                        ['Power Required for Compressor is ',round(Compressorpower/1000,4),'KW'],
                        ['Mass Flow Rate(input paramerter)of VCR',mdot,'Kg/s']]
                print(tabulate(table2,headers=['Parameter','Value','Unit']))
                print('--------------------------------End of Vapor Compression Cycle--------------------------------')
                print('\n')
                print('--------------------------------Organic Rankine Cycle--------------------------------')
                table3=[['a',Ta,round(Sa/1000,4),round(Ha/1000,4),Pa/1e5],
                        ['a_',Ta,round(S_a/1000,4),round(H_a/1000,4),Pa/1e5],
                        ['b',Tb,round(Sb/1000,4),round(Hb/1000,4),Pb/1e5],
                        ['c',Tc,round(Sc/1000,4),round(Hc/1000,4),Pc/1e5],
                        ['d',Td,round(Sd/1000,4),round(Hd/1000,4),Pd/1e5],
                        ['e',Te,round(Se/1000,4),round(He/1000,4),Pe/1e5],
                        ['e_',Te,round(S_e/1000,4),round(H_e/1000,4),Pe/1e5],
                        ['f',Tf,round(Sf/1000,4),round(Hf/1000,4),Pf/1e5],
                        ['g',Tg,round(Sg/1000,4),round(Hg/1000,4),Pg/1e5]]
                print(tabulate(table3,headers=['State Point','Temperature(in K)','Entropy(S)(in KJ/Kg)','Enthaply(H)(KJ/Kg)','Pressure(bar)']))
                table4=[['Qin from VCR ',round(Qin/1000,4),'KJ/Kg'],
                        ['Turbine Work Ideal',round(Wt/1000,4),'KW'],
                        ['Pump Work Ideal',round(Wp/1000,4),'KW'],
                        ['ORC Efficiency',ORC_Efficiency,'%'],
                        ['Actual Turbine Work ',round(Wtactual/1000,4),'KW'],
                        ['Actual Pump Work',round(Wpactual/1000,4),'KW'],
                        ['Actual ORC Efficiency is',ORC_Efficiency_actual,'%'],
                        ['Mass Flow Rate of Working Fluid in ORC',mdot_ORC,'Kg/s'],
                        ['Specific Workinf Fluid Consumption',SteamRate,'Kg/KWH']]
                print(tabulate(table4, headers=['Parameter', 'Value', 'Unit']))
                data=pd.read_excel('./EXCELTEST-1.xlsx',index_col=[0])
                finaldata=data.append({'Working Fluid in VCR':fluid,
                                'Mass Flow Rate(input in Kg/s)':mdot,
                                'Evaporator Temperature(C)':T1-273.15,
                                'SuperHeating Temperature(C)':T2-273.15,
                                'Condensor Temperature(C)':T4-273.15,
                                'SubCooling Temperature(C)':T6-273.15, 
                                'COPIdeal':COPideal,
                                'COPActual':COPactual,
                                'Refrigeration Capacity(tons)':RC,
                                'QoutcondensorVCR(KW)':round(CondesorHeat/1000,4),
                                'CompressorPower(KW)':round(Compressorpower/1000,4),
                                'QinfromVCR(KW)':round(Qin/1000),
                                'Working Fluid in ORC':fluid2,
                                'Mass Flow Rate of Working Fluid in ORC':mdot_ORC,
                                'Pressure of Boiler working Fluid(ORC)(bar)':Pa/1e5,
                                'Condensor Pressure of Working Fluid(ORC)(bar)':Pe/1e5,
                                'Turbine Work Ideal(KW)':round(Wt/1000,4),
                                'Pump Work Ideal(KW)':round(Wp/1000,4),
                                'ORC Efficiency(%)':ORC_Efficiency,
                                'Actual Turbine Work(KW)':round(Wtactual/1000,4),
                                'Actual Pump Work(KW)':round(Wpactual/1000,4),
                                'Actual ORC Efficiency is(%)':ORC_Efficiency_actual,
                                'Specific Workinf Fluid Consumption(Kg/KWH)':SteamRate},ignore_index=True)
                finaldata=pd.DataFrame(finaldata)
                writer=pd.ExcelWriter('EXCELTEST-1.xlsx')
                finaldata.to_excel(writer)
                writer.save()
                print('#####################################################')
                print('Data Saved in EXCEL')
                print('#####################################################')
                #VCR Plotting
                Tcritical_VCR=PropsSI("Tcrit",fluid)
                T_vcr = np.arange(180, Tcritical_VCR,0.1)
                Sf_vcr = PropsSI('S','T',T_vcr,'Q',0,fluid)/1000
                Sg_vcr= PropsSI('S','T',T_vcr,'Q',1,fluid)/1000
                plt.subplot(2,2,1)
                plt.title("TS Diagram")
                plt.xlabel("X axis(Entropy in KJ/Kg)")
                plt.ylabel("Y axis(Temperature in K)")
                plt.plot(Sg_vcr,T_vcr,color='black')
                plt.plot(Sf_vcr,T_vcr,color='black')
                plt.annotate("1",(round(S1/1000,4),T1),fontsize=16)
                plt.annotate("2",(round(S2/1000,4),T2),fontsize=16)
                plt.annotate("3",(round(S3/1000,4),T3),fontsize=16)
                plt.annotate("4",(round(S4/1000,4),T4),fontsize=16)
                plt.annotate("5",(round(S5/1000,4),T5),fontsize=16)
                plt.annotate("6",(round(S6/1000,4),T6),fontsize=16)
                plt.annotate("7",(round(S7/1000,4),T7),fontsize=16)
                plt.scatter([round(S1/1000,2),round(S2/1000,2),round(S3/1000,2),round(S4/1000,2),round(S5/1000,2),round(S6/1000,2),round(S7/1000,2)],[T1,T2,T3,T4,T5,T6,T7],color="b")
                plt.plot([round(S1/1000,2),round(S2/1000,2),round(S3/1000,2),round(S4/1000,2),round(S5/1000,2),round(S6/1000,2),round(S7/1000,2),round(S1/1000,2)],[T1,T2,T3,T4,T5,T6,T7,T1],color='red')
                Pcritical_VCR=PropsSI("Pcrit",fluid)
                P_vcr= np.arange(.1e5, Pcritical_VCR,0.235e5)
                Hf_vcr = PropsSI('H','P',P_vcr,'Q',0,fluid)/1000
                Hg_vcr= PropsSI('H','P',P_vcr,'Q',1,fluid)/1000
                P_vcr=P_vcr/100000
                plt.subplot(2,2,2)
                plt.title("PH Diagram")
                plt.xlabel("X axis(Enthaply in KJ/Kg)")
                plt.ylabel("Y axis(Pressure in Bar)")
                plt.plot(Hf_vcr,P_vcr,color='black')
                plt.plot(Hg_vcr,P_vcr,color='black')
                plt.annotate("1",(round(H1/1000,4),P1),fontsize=16)
                plt.annotate("2",(round(H2/1000,4),P2),fontsize=16)
                plt.annotate("3",(round(H3/1000,4),P3),fontsize=16)
                plt.annotate("4",(round(H4/1000,4),P4),fontsize=16)
                plt.annotate("5",(round(H5/1000,4),P5),fontsize=16)
                plt.annotate("6",(round(H6/1000,4),P6),fontsize=16)
                plt.annotate("7",(round(H7/1000,4),P7),fontsize=16)
                plt.scatter([round(H1/1000,2),round(H2/1000,2),round(H3/1000,2),round(H4/1000,2),round(H5/1000,2),round(H6/1000,2),round(H7/1000,2)],[P1,P2,P3,P4,P5,P6,P7],color="b")
                plt.plot([round(H1/1000,2),round(H2/1000,2),round(H3/1000,2),round(H4/1000,2),round(H5/1000,2),round(H6/1000,2),round(H7/1000,2),round(H1/1000,2)],[P1,P2,P3,P4,P5,P6,P7,P1],color='red')
                #ORC PLotting
                Tcritical_ORC=PropsSI("Tcrit",fluid2)
                T_ORC=np.arange(170, Tcritical_ORC,0.1)
                Sf_ORC = PropsSI('S','T',T_ORC,'Q',0,fluid2)/1000
                Sg_ORC= PropsSI('S','T',T_ORC,'Q',1,fluid2)/1000
                plt.subplot(2,2,3)
                plt.title("TS Diagram")
                plt.xlabel("X axis(Entropy in KJ/Kg)")
                plt.ylabel("Y axis(Temperature in K)")
                plt.plot(Sg_ORC,T_ORC,color='black')
                plt.plot(Sf_ORC,T_ORC,color='black')
                plt.annotate("a",(round(Sa/1000,4),Ta),fontsize=16)
                plt.annotate("a'",(round(S_a/1000,4),Ta),fontsize=16)
                plt.annotate("b",(round(Sb/1000,4),Tb),fontsize=16)
                plt.annotate("c",(round(Sc/1000,4),Tc),fontsize=16)
                plt.annotate("d",(round(Sd/1000,4),Td),fontsize=16)
                plt.annotate("e",(round(Se/1000,4),Te),fontsize=16)
                plt.annotate("e'",(round(S_e/1000,4),Te),fontsize=16)
                plt.annotate("f",(round(Sf/1000,4),Tf),fontsize=16)
                plt.annotate("g",(round(Sg/1000,4),Tg),fontsize=16)
                plt.scatter([round(Sa/1000,2),round(S_a/1000,2),round(Sb/1000,2),round(Sc/1000,2),round(Sd/1000,2),round(Se/1000,2),round(S_e/1000,2),round(Sf/1000,2),round(Sg/1000,2)],[Ta,Ta,Tb,Tc,Td,Te,Te,Tf,Tg],color="b")
                plt.plot([round(Sa/1000,2),round(Sb/1000,2),round(Sc/1000,2),round(Sd/1000,2),round(Se/1000,2),round(Sf/1000,2),round(Sg/1000,2),round(Sa/1000,2)],[Ta,Tb,Tc,Td,Te,Tf,Tg,Ta],color='red')
                plt.plot([round(Sd/1000,2),round(S_e/1000,2),round(Se/1000,2)],[Td,Te,Te],'--')
                plt.plot([round(Sg/1000,2),round(S_a/1000,2),round(Sb/1000,2)],[Tg,Ta,Tb],'--',color='blue')
                Pcritical_ORC=PropsSI("Pcrit",fluid)
                P_ORC= np.arange(.1e5,Pcritical_ORC,0.285e5)
                Hf_ORC = PropsSI('H','P',P_ORC,'Q',0,fluid2)/1000
                Hg_ORC= PropsSI('H','P',P_ORC,'Q',1,fluid2)/1000
                P_ORC=P_ORC/1e5
                plt.subplot(2,2,4)
                plt.title("PH Diagram")
                plt.xlabel("X axis(Enthaply in KJ/Kg)")
                plt.ylabel("Y axis(Pressure in Bar)")
                plt.plot(Hf_ORC,P_ORC,color='black')
                plt.plot(Hg_ORC,P_ORC,color='black')
                plt.annotate("a",(round(Ha/1000,4),Pa/1e5),fontsize=16)
                plt.annotate("a'",(round(H_a/1000,4),Pa/1e5),fontsize=16)
                plt.annotate("b",(round(Hb/1000,4),Pb/1e5),fontsize=16)
                plt.annotate("c",(round(Hc/1000,4),Pc/1e5),fontsize=16)
                plt.annotate("d",(round(Hd/1000,4),Pd/1e5),fontsize=16)
                plt.annotate("e",(round(He/1000,4),Pe/1e5),fontsize=16)
                plt.annotate("e'",(round(H_e/1000,4),Pe/1e5),fontsize=16)
                plt.annotate("f",(round(Hf/1000,4),Pf/1e5),fontsize=16)
                plt.annotate("g",(round(Hg/1000,4),Pg/1e5),fontsize=16)
                plt.scatter([round(Ha/1000,2),round(H_a/1000,2),round(Hb/1000,2),round(Hc/1000,2),round(Hd/1000,2),round(He/1000,2),round(H_e/1000,2),round(Hf/1000,2),round(Hg/1000,2)],[Pa/1e5,Pa/1e5,Pb/1e5,Pc/1e5,Pd/1e5,Pe/1e5,Pe/1e5,Pf/1e5,Pg/1e5],color="b")
                plt.plot([round(Ha/1000,2),round(Hb/1000,2),round(Hc/1000,2),round(Hd/1000,2),round(He/1000,2),round(Hf/1000,2),round(Hg/1000,2),round(Ha/1000,2)],[Pa/1e5,Pb/1e5,Pc/1e5,Pd/1e5,Pe/1e5,Pf/1e5,Pf/1e5,Pa/1e5],color='red')
                plt.plot([round(Hd/1000,2),round(H_e/1000,2),round(He/1000,2)],[Pd/1e5,Pe/1e5,Pe/1e5],'--')
                plt.plot([round(Hg/1000,2),round(H_a/1000,2),round(Hb/1000,2)],[Pg/1e5,Pa/1e5,Pb/1e5],'--',color='blue')
                plt.show()
                input('Press Enter Key to Exit')     
        if x>1:
            print('Exit Quality of Fluid from Turbine is SuperHeated (i.e X>1)')
            Te=PropsSI('T','P',Pg,'S',Se,fluid2)
            He=round(PropsSI('H','P',Pg,'T',Te,fluid2),4)
            #He=round(PropsSI('H','P',Hg,'S',Se,fluid2),4)
            Hf=round(PropsSI('H','P',Pf,'Q',1,fluid2),4)
            Sf=round(PropsSI('S','P',Pf,'Q',1,fluid2),4)
            condensorORC(fluid,fluid2,Pa,Pb,Pc,Pd,Pe,Pf,Pg,Ha,Hb,Hc,Hd,He,Hf,Hg,Tb,Tc,Td,Te,Tf,Tg,Sb,Sc,Sd,Se,Sf,T1,T2,T3,T4,T5,T6,T7,S1,S2,S3,S4,S5,S6,S7,P1,P2,P3,P4,P5,P6,P7,H1,H2,H3,H4,H5,H6,H7,COPideal,COPactual,RC,CondesorHeat,Compressorpower,mdot)
        else:
            Sf=Sd
            Te=Tf
            Hg1=round(PropsSI('H','P',Pg,'Q',1,fluid2),4)
            He=Hf=(Hg+x*(Hg1-Hg))
            condensorORC(fluid,fluid2,Pa,Pb,Pc,Pd,Pe,Pf,Pg,Ha,Hb,Hc,Hd,He,Hf,Hg,Tb,Tc,Td,Te,Tf,Tg,Sb,Sc,Sd,Se,Sf,T1,T2,T3,T4,T5,T6,T7,S1,S2,S3,S4,S5,S6,S7,P1,P2,P3,P4,P5,P6,P7,H1,H2,H3,H4,H5,H6,H7,COPideal,COPactual,RC,CondesorHeat,Compressorpower,mdot)
        
    if Hd<Hc:
        print('Quality of Steam at Inlet to Turbine is not dry (x<1)')
        print('Therefore it is compensated by reducing mass flow flow rate from intially assumed 1Kg/s and required mass flow rate(ORC) is calulated in order to reach at leat saturated state')
        print('Do you want to Proceed by Compensating the mass flow rate of ORC working Fluid?')
        opt=int(input('If Yes.Plz Enter 1\n'))
        if opt==1:
            mdot_ORC=round(Qin/(Hc-Ha),3)
            Hd=Hc
            Sd=Sc
            Se=Sd
            Td=Tc
            ORC_sat(Tc,Td,Hc,Hd,Sc,Sd,Se,Sf0,Sf1,mdot_ORC)
        else:
            print('You choosen to Termite the Program')
    else:
        Cc=round(PropsSI('C','P',Pc,'Q',1,fluid2),4)
        #Tsup=((Hd-Hc)/Cc)+Tc
        Tsup=round(PropsSI('T','P',Pd,'H',Hd,fluid2),4)
        Td=Tsup
        Sd=round(PropsSI('S','P',Pd,'T',Td,fluid2),4)
        Se=Sd
        mdot_ORC=1
        ORC_sat(Tc,Td,Hc,Hd,Sc,Sd,Se,Sf0,Sf1,mdot_ORC) 
print('-------------VAPOUR COMPRESSION REFRIGERATION CYCLE------------------')
print('Select any one of the Option')
option_a=int(input('1)No SuperHeating and No SubCooling\n'
      '2)Only SuperHeating No SubCooling\n'
      '3)Only Subcooling No SuperHeating\n'
      '4)Both SuperHeating and SubCooling\n'))
          
if option_a==1:
    fluid=input('Enter the Valid Refrigerent Name\n')
    print("#########")
    print("NOTE:",PropsSI("Tcrit",fluid)-273,"is the maximum Temperature Allowed")
    print("#########")
    mdot=float(input('Enter the mass flow rate in Kg/s\n'))
    T1=float(input('Enter the Evaporator Temperature in Celsius\n'))
    T4=float(input('Enter the COndensor temperature in Celius\n'))
    T2=T1
    T6=T4

    work(fluid,mdot,T1,T4,T2,T6)
    
elif option_a==2:
    fluid=input('Enter the Valid Refrigerent Name\n')
    print("#########")
    print("NOTE:",PropsSI("Tcrit",fluid)-273,"is the maximum Temperature Allowed")
    print("#########")
    mdot=float(input('Enter the mass flow rate in Kg/s\n'))
    T1=float(input('Enter the Evaporator Temperature in Celsius\n'))
    T4=float(input('Enter the COndensor temperature in Celius\n'))
    T2=float(input('Enter the Superheat temperature in Celius\n'))
    T6=T4

    work(fluid,mdot,T1,T4,T2,T6)
   
elif option_a==3:
    fluid=input('Enter the Valid Refrigerent Name\n')
    print("#########")
    print("NOTE:",PropsSI("Tcrit",fluid)-273,"is the maximum Temperature Allowed")
    print("#########")
    mdot=float(input('Enter the mass flow rate in Kg/s\n'))
    T1=float(input('Enter the Evaporator Temperature in Celsius\n'))
    T4=float(input('Enter the COndensor temperature in Celius\n'))
    T2=T1
    T6=float(input('Enter the SubCooling temperature in Celius\n'))
    
    work(fluid,mdot,T1,T4,T2,T6)
 
elif option_a==4:
    fluid=input('Enter the Valid Refrigerent Name\n')
    print("#########")
    print("NOTE:",PropsSI("Tcrit",fluid)-273,"is the maximum Temperature Allowed")
    print("#########")
    mdot=float(input('Enter the mass flow rate in Kg/s\n'))
    T1=float(input('Enter the Evaporator Temperature in Celsius\n'))
    T4=float(input('Enter the COndensor temperature in Celius\n'))
    T2=float(input('Enter the Superheat temperature in Celius\n'))
    T6=float(input('Enter the SubCooling temperature in Celius\n'))

    work(fluid,mdot,T1,T4,T2,T6)

else:    
    print('Error! Plz Try again.')



            
    

