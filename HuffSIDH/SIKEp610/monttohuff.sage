pp=(2^372)*(3^239)-1
F=GF(pp)
K.<i>=GF(pp^2,modulus=x^2+1)


# montgomery point for Alice
#xQ20+xQ21*i
xQ20=2141306817115640038993862181765492931500972641785085350947220052387904285572880641340973064166142567836013875254694387171253686309941716659316458912387577733949836538957913038827049600784199475852241950366745525218545897974030
xQ21=3462139332946457351613256472394206428708902796442953775366163334092154110965882738568198703970873585433373209361346635159539572450380354861849053160099274714477389537230240941950694098078938506114926490884815238642940114646471

xAQ=xQ20+xQ21*i

t0=K(1/xAQ)

print('Alice Q on Huff coordinate=',t0)

#xP20+xP21*i

xP20=6392653750389499144594798145926583444354409039583606858359849849952175922439864668710005278778374444865434395866810086432449384862893379044515952989369610078853427396744868696171536543536691791447062025853410930611265727936250
xP21=1993436594774325849134512059294984900092784525889309731865924302516333954639088053874784298502590615500299169775920492705338536494948131816798689951115100682529755597000522266521373993031102110056636212011317531431383151419457

xAP=xP20+xP21*i


t0=K(1/xAP)

print('Alice P on Huff coordinate=',t0)


#xR20+xR21
xR20=8920767099889265354270903767976190589347780532338110811601966065539889197224092469005030134458217535515117192741315988190233956060833229415177951494574325212388096836354324553840374916089004540013556866495215628182837045568955
xR21=7485054888408186199988187951947272686736588976237807119427933559723580777349951819085282578750816796140208081851358534302244661919771165437390420226722284298495663460285754671947773474537816261384790115136982920513335294841189

xAR=xR20+xR21*i


t0=K(1/xAR)

print('Alice R on Huff coordinate=',t0)



#montgomery point for bob
#xQ30+xQ31*i

xQ30=8511019311184218364369646557252415173765495706553816075924705665118452233626065646710473783622713033004291163624815196837285573011460225202480327659675906178097793861901498715713456623486950717478378070640744966049319812387057

xBQ=xQ30


t0=K(1/xBQ)

print('Bob Q on Huff coordinate=',t0)


#xP30
xP30=8917296521312037934018779383955573993264347381711746006893090701936594237023232166684104618705889709151882479004705477580016078091729035680421189641116128972780575456542140523559278016360828837041023736833315936592076152779518
     
xBP=xP30

t0=K(1/xBP)

print('Bob P on Huff coordinate=',t0)


#xR30+xR31*i
xR30=7948472920477987382026230091597296324770743133322382953559027144948810258767547973629670019619498810386907603969957459770220385943384654803338204542090155222002874700542244358439397441360179573140328384790672595434635205779427
xR31=8399745264634857705578681279266694258039523836471824133317786739845301013095107323108209992417182837952421722872170331095872127782294209276940273704388772866734580589843768351599575913837321305038643830671122310284823719495713

xBR=xR30+xR31*i

t0=K(1/xBR)

print('Bob R on Huff coordinate=',t0)




## To mont form

AQ1=2160705014733594267295946005265155614468510711066162509895591195578879937601239815053181288152754084500387263246189419374891161519008296549071364329798210535377727017843628054763566387822825387902093936433450398690487281528308
AQ0=1880071532701365469674459205645396243871993547553265366888595082145479048767590861267916676590930922208016842640939042989076219889149840428108113167954127269825431695693849319269077320862629689293980867610914396903113292926297

AP1=4642458033700095254072871244277677748214783002041880889458919399546656795415037394875739677188725472232771185837888192756303618933512776894579821414919433715652000231913161316882306955379351895004616904860622979167314650956839
AP0=5855777938868749393157097443941145331520565631621986746061635845093895539088668766698704356537301802159014031737157459903980848772525023885124537747501745071282225033848836733133511109316133997291204176450625227235456934025523

AR1=1432024966259981857940126279509966478350565160171052118797518405622645922375450695881971940608728061412179784758421153726074873172675629748363440347051173619646628748365425512529290011749023695302363024157664701584107131486018
AR0=364709152190734789183543668560098560247527500818347199713623604613747083213335299410455825903650447355871003634234665013077491286891359685755732814899256141893748556030442410327789650357699345209846138400374226404851544159103


BQ0=3730898147656993489396481267626659156206944759219145400675841807903397874410894185421783064234663096267843054041210031419441525568857379526531733919764486296275273416305528293538997019543433371542809649946485704814427902230142


BP0=8138763166322702834602773966715226729460438343085860906768535726848297161622348191824873988691958335658534890594753009007054122428530563388297458376742375434624070755962310536072020734698278179063689766461804956400805622081295

    

BR1=7722539722710165823804275554317341027883689907040189755045605640246739141146203882891924702460659327183310776970648399387302706037332860029677819799635814035192964621255516829602088314196746595670459802572192632928303646718460
BR0=5536512464204382375865442413971438070515754169282094321226060150155657471477342741371504541181042412024608996252824634690188662126080875811441259015951052708793574592998779859135118175743238049423215181377905594008649873376537


aq1=mod(AQ1*(2^768),pp)
aq0=mod(AQ0*(2^768),pp)

ap1=mod(AP1*(2^768),pp)
ap0=mod(AP0*(2^768),pp)

ar1=mod(AR1*(2^768),pp)
ar0=mod(AR0*(2^768),pp)

bq0=mod(BQ0*(2^768),pp)

bp0=mod(BP0*(2^768),pp)

br1=mod(BR1*(2^768),pp)
br0=mod(BR0*(2^768),pp)

print('')
print('')
print('')

print('Alice P0 on Huff coordinate=',hex(ap0))
print('Alice P1 on Huff coordinate=',hex(ap1))
print('Alice Q0 on Huff coordinate=',hex(aq0))
print('Alice Q1 on Huff coordinate=',hex(aq1))
print('Alice R0 on Huff coordinate=',hex(ar0))
print('Alice R1 on Huff coordinate=',hex(ar1))

print('Bob BQ0 on Huff coordinate=',hex(bq0))
print('Bob BP0 on Huff coordinate=',hex(bp0))
print('Bob BR0 on Huff coordinate=',hex(br0))
print('Bob BR1 on Huff coordinate=',hex(br1))