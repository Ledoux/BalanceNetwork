#/##################/#
# Architecture
#

#NeuronGroups
NeuronGroupStrsList=map(
    lambda __NeuronTypeStr:
    __NeuronTypeStr+'NeuronGroup',
    NeuronTypeStrsList
)

#Inputs
PoissonInputStrsList=map(
    lambda __NeuronTypeStr:
    __NeuronTypeStr+'PoissonInput',
    NeuronTypeStrsList
)

#Synapses
import itertools
NeuronTypeStrsTuplesList=list(
    itertools.product(
        NeuronTypeStrsList,
        NeuronTypeStrsList
    )
)
SynapsesStrsList=map(
    lambda __NeuronTypeStrsTuple:
    __NeuronTypeStrsTuple[0]+'To'+__NeuronTypeStrsTuple[1]+'Synapses',
    NeuronTypeStrsTuplesList
)

#StateMonitor
StateMonitorStrsList=map(
    lambda __NeuronTypeStr:
    __NeuronTypeStr+'StateMonitor',
    NeuronTypeStrsList
)

#SpikeMonitor
SpikeMonitorStrsList=map(
    lambda __NeuronTypeStr:
    __NeuronTypeStr+'SpikeMonitor',
    NeuronTypeStrsList
)

#PopulationRateMonitor
PopulationRateMonitorStrsList=map(
    lambda __NeuronTypeStr:
    __NeuronTypeStr+'PopulationRateMonitor',
    NeuronTypeStrsList
)

#/##################/#
# Brian import
#

#import
import brian2

#derive
class SynapsesClass(brian2.Synapses):
    
    def connect(self,*_LiargVariablesList,**_KwargVariablesDict):
        brian2.Synapses.connect(self,*_LiargVariablesList,**_KwargVariablesDict)
        return self
    
    def setWeight(self,_WeigthFloatsArray):
        self.J[:]=_WeigthFloatsArray.reshape(
                self.source.N*self.target.N
            )
        return self

#/##################/#
# Reset Network and all brian objects
#

#build network
InstanceTrackerSet=brian2.Nameable.__instances__()
map(
    lambda __BrianVariable:
    InstanceTrackerSet.remove(__BrianVariable),
    list(InstanceTrackerSet)
);
    
#build network
BalanceNetwork=brian2.Network()

#/##################/#
# build neuron groups
#

#map
map(
    lambda __NeuronGroupStr:
    BalanceNetwork.add(
        brian2.NeuronGroup(
            N=UnitsInt
            if __NeuronGroupStr.startswith('Exc')
            else InhRatioFloat*UnitsInt,
            model='''
                mu : volt
                dv/dt= ( -(v-('''+str(RestFloat)+'''*mV)) + mu + xi*sqrt('''+str(ConstantTimeFloat
                )+'''*ms)*0.*mV ) /('''+str(ConstantTimeFloat)+'''*ms) : volt  (unless refractory)
            ''',
            threshold="v>"+str(ThresholdFloat)+"*mV",
            reset="v="+str(ResetFloat)+"*mV",
            refractory=str(RefractoryTimeFloat)+"*ms",
            name=__NeuronGroupStr
        )
    ) 
    if __NeuronGroupStr not in BalanceNetwork
    else None,
    NeuronGroupStrsList
);

#set
BalanceNetwork['ExcNeuronGroup'].clock.dt=RunStepTimeFloat*brian2.ms

#build poisson inputs
map(
    lambda __PoissonInputStr,__NeuronGroupStr:
    BalanceNetwork.add(
        brian2.PoissonInput(
            BalanceNetwork[__NeuronGroupStr],
            'v',
            SparsityFloat*UnitsInt,
            PoissonFrequencyFloat*brian2.Hz,
            weight=ExcWeightFloat*brian2.mV
        )
    ),
    PoissonInputStrsList,
    NeuronGroupStrsList
);

#/##################/#
# build synapses
#

#import
import scipy.stats

#Check
if LateralIsBool:

    #map
    map(
        lambda __NeuronTypeStrsTuple,__SynapsesStr:
        BalanceNetwork.add(
            SynapsesClass(
                BalanceNetwork[__NeuronTypeStrsTuple[0]+'NeuronGroup'],
                BalanceNetwork[__NeuronTypeStrsTuple[1]+'NeuronGroup'],
                model='''
                    J:1
                ''',
                pre="v+=J*mV",
                delay=DelayTimeFloat*brian2.ms,
                name=__SynapsesStr
            ).connect(
                LateralIsBool
            ).setWeight(
                ExcWeightFloat*scipy.stats.bernoulli.rvs(
                    SparsityFloat,
                    size=(
                        BalanceNetwork[__NeuronTypeStrsTuple[0]+'NeuronGroup'].N,
                        BalanceNetwork[__NeuronTypeStrsTuple[1]+'NeuronGroup'].N
                    )
                )
                if __NeuronTypeStrsTuple[0]=='Exc'
                else
                -InhScaleFloat*ExcWeightFloat*scipy.stats.bernoulli.rvs(
                    SparsityFloat,
                    size=(
                        BalanceNetwork[__NeuronTypeStrsTuple[0]+'NeuronGroup'].N,
                        BalanceNetwork[__NeuronTypeStrsTuple[1]+'NeuronGroup'].N
                    )
                )
            )
        )
        if __SynapsesStr not in BalanceNetwork
        else None,
        NeuronTypeStrsTuplesList,
        SynapsesStrsList
    );

#/##################/#
# build Monitors
#

#build state monitor for v
map(
    lambda __NeuronGroupStr,__StateMonitorStr:
    BalanceNetwork.add(
        brian2.StateMonitor(
            BalanceNetwork[__NeuronGroupStr],
            'v',
            [0,1],
            name=__StateMonitorStr,
        )
    )
    if __StateMonitorStr not in BalanceNetwork
    else None,
    NeuronGroupStrsList,
    StateMonitorStrsList
);

#build spike monitor
map(
    lambda __NeuronGroupStr,__SpikeMonitorStr:
    BalanceNetwork.add(
        brian2.SpikeMonitor(
            BalanceNetwork[__NeuronGroupStr],
            name=__SpikeMonitorStr
        )
    )
    if __SpikeMonitorStr not in BalanceNetwork
    else None,
    NeuronGroupStrsList,
    SpikeMonitorStrsList
);

#build rate monitor
map(
    lambda __NeuronGroupStr,__PopulationRateMonitorStr:
    BalanceNetwork.add(
        brian2.PopulationRateMonitor(
            BalanceNetwork[__NeuronGroupStr],
            name=__PopulationRateMonitorStr
        )
    )
    if __PopulationRateMonitorStr not in BalanceNetwork
    else None,
    NeuronGroupStrsList,
    PopulationRateMonitorStrsList
);

#/##################/#
# Random initial values 
#

#map
map(
    lambda __NeuronGroupStr:
    setattr(
        BalanceNetwork[__NeuronGroupStr],
        'v',
        (
            RestFloat+scipy.stats.uniform.rvs(
                size=BalanceNetwork[__NeuronGroupStr].N
            )
        )*brian2.mV
    ),
    NeuronGroupStrsList
)

#/##################/#
# Run 
#

#run
BalanceNetwork.run(RunDurationTimeFloat*brian2.ms);

#/##################/#
# Plot
#

#dict
ColorDict={
	"Exc":"blue",
	"Inh":"red"
}

#import 
from matplotlib import pyplot

#init
RunFigure=pyplot.figure(
    figsize=(
            15,10
    )
)


#/##################/#
# Voltage Traces
#

#init
StateAxes=RunFigure.add_subplot(311) 

#set
XlimList=[
        BalanceNetwork['ExcStateMonitor'].t[0],
        BalanceNetwork['ExcStateMonitor'].t[-1]
]

#map
map(
    lambda __NeuronTypeStr,__StateMonitorStr:
    StateAxes.plot(
            BalanceNetwork[__StateMonitorStr].t,
            BalanceNetwork[__StateMonitorStr].v.T,
            color=ColorDict[__NeuronTypeStr],
            linewidth=3
        ),
    NeuronTypeStrsList,
    StateMonitorStrsList
)
StateAxes.plot(XlimList,[0.001*ThresholdFloat]*2,'--',color='black')
StateAxes.plot(XlimList,[0.001*RestFloat]*2,'--',color='black')
StateAxes.plot(XlimList,[0.001*ResetFloat]*2,'--',color='black')

#/##################/#
# Spike Traces
#

#init
SpikeAxes=RunFigure.add_subplot(312) 

#map
map(
    lambda __IndexInt,__NeuronTypeStr,__SpikeMonitorStr:
    SpikeAxes.plot(
            BalanceNetwork[__SpikeMonitorStr].t,
            BalanceNetwork[__SpikeMonitorStr].i+(
                BalanceNetwork[
                    NeuronTypeStrsList[
                        __IndexInt-1
                    ]+'NeuronGroup'
                ].N
                if __IndexInt>0 else 0.
            ),
            linestyle='',
            marker='o',
            markersize=2,
            color=ColorDict[__NeuronTypeStr],
            label='$'+__NeuronTypeStr+'$'
        ),
    xrange(len(NeuronTypeStrsList)),
    NeuronTypeStrsList,
    SpikeMonitorStrsList
)
SpikeAxes.set_xlim(XlimList)
SpikeAxes.set_ylim([0,UnitsInt*(1.+InhRatioFloat)])
SpikeAxes.legend()


#/##################/#
# Rate Traces
#

#init
RateAxes=RunFigure.add_subplot(313) 

#set bins
WindowTimeFloat = 0.4
WindowLengthInt = int(WindowTimeFloat*brian2.ms/brian2.defaultclock.dt)
CumsumArraysList=map(
    lambda __PopulationRateMonitorStr:
    numpy.cumsum(numpy.insert(BalanceNetwork[__PopulationRateMonitorStr].rate,0,0)),
    PopulationRateMonitorStrsList
)
BinRateArraysList=map(
    lambda __CumsumArray:
    (__CumsumArray[WindowLengthInt:] - __CumsumArray[:-WindowLengthInt]) / WindowLengthInt,
    CumsumArraysList
)

#map
map(
    lambda __NeuronTypeStr,__PopulationRateMonitorStr,__BinRateArray:
    RateAxes.plot(
            BalanceNetwork[__PopulationRateMonitorStr].t[WindowLengthInt-1:]/brian2.ms,
            __BinRateArray,
            color=ColorDict[__NeuronTypeStr],
            linewidth=3
    ),
    NeuronTypeStrsList,
    PopulationRateMonitorStrsList,
    BinRateArraysList
)


