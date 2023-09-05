use std::str::FromStr;

use crate::FieldGC;
use num_bigint::BigUint;

use crate::PoseidonConstants;

// We dynamically set the constants for the secp256k1 curve instead a hardcoding for convenience.
// Hardcoding requires us to use the `ark_ff::Fq` type, which
// is hard to use in combination with other generic types.

pub fn secp256k1<F: FieldGC>() -> PoseidonConstants<F> {
    let num_full_rounds = 8;
    let num_partial_rounds = 56;

    let mds_matrix = vec![
        [
            "92469348809186613947252340883344274339611751744959319352506666082431267346705",
            "100938028378191533449096235266991198229563815869344032449592738345766724371160",
            "77486311749148948616988559783475694076613010381924638436641318334458515006661",
        ]
        .map(|y| F::from(BigUint::from_str(y).unwrap()))
        .to_vec(),
        [
            "110352262556914082363749654180080464794716701228558638957603951672835474954408",
            "27607004873684391669404739690441550149894883072418944161048725383958774443141",
            "29671705769502357195586268679831947082918094959101307962374709600277676341325",
        ]
        .map(|y| F::from(BigUint::from_str(y).unwrap()))
        .to_vec(),
        [
            "77762103796341032609398578911486222569419103128091016773380377798879650228751",
            "1753012011204964731088925227042671869111026487299375073665493007998674391999",
            "70274477372358662369456035572054501601454406272695978931839980644925236550307",
        ]
        .map(|y| F::from(BigUint::from_str(y).unwrap()))
        .to_vec(),
    ];

    let round_keys = [
        "15180568604901803243989155929934437997245952775071395385994322939386074967328",
        "98155933184944822056372510812105826951789406432246960633912199752807271851218",
        "32585497418154084368870158853355239726261349829448673320273043226636389078017",
        "66713968576806622579829258440960693099797917756640662361943757758980796487698",
        "61296025743283504825054745787375839406507895949474930140819919915792438454216",
        "64548089412749542282115556935384382035671782881737715696939837764375912217104",
        "108421562972909537718478936575770973463273651828765393113349044862621092658552",
        "93957623861448681916560847065407918286434708744548934125771289238599801659600",
        "31886767595881910145119755249133120645312710313371225820300496900248094187131",
        "36511615103248888903406040506250394762206798360602726106046630438239169384653",
        "21193239787133737740669439860809806837993750509086389566475677877580362491125",
        "15159189447883181997488877417695825734356570617827322308691834229181804753656",
        "19272373877630561389686073945290625876718814210798194797601715657476609730306",
        "23132197996397121955527964729507651432518694856862854469217474256539272053037",
        "9869753235007825662020275771343858285582964429845049469800863115040150206544",
        "36536341316285671890133896506951910369952562161551585116256678375995315827743",
        "62582239167707347777855528698896708360409296899261565735324151945083720570858",
        "96597358901965097853721114962031771931271685249979807653919643952343419105640",
        "99475971754252188104003224702005940217163363685728394033034788135108600073953",
        "52080483875928847502018688921126796935417602445765802481027972679966274137987",
        "101922748752417217354391348649359865075718358385248454632698502400961567227929",
        "26980595292132221181330746499613907829041623688147011560382352796984836870749",
        "7059991836806083192408106370472821784612460308866802565871813230060135266390",
        "19329812920723038526370491239817117039289784665617181727933894076969997926129",
        "65570620823578601926240439251563587376966657231502120214692324496443514623818",
        "58403733332589349613112270854204921427257113546270812628317365115158685715742",
        "45021021211732634759643776743541935700591354899980928498981462362035961745443",
        "313468157086800401026946312285365733155132234906935411743639256319782592571",
        "101316949793045093761117346380310841944294663456931203380573537653884068660109",
        "23683935571424619534194393788101669168630123784066421490798386323411538828592",
        "45470730427236677197026094498490008082250264942279323465121581539984407294442",
        "48141067373531800337373447278127981363951468257064369512416205750641258258193",
        "42554919225040466028330117313396362347164995917041931400909482795914116747618",
        "11551941832988244108260444347046942236051939264069344013774353630451796870907",
        "60185799182545404739011626517355854847787627814101363386450657535504094743765",
        "81823160578900678880708744457872721685515019032370491632046212317701226128393",
        "7165646831054215773988859638722974820791178194871546344315162343128362695647",
        "75289707601640398424243937567716657896680380639974371761136292031415717685949",
        "7150842764562742184396161198129263121409208675362553300851082062734889620953",
        "24380904705269761063866540342138412601132455197711667167747524315310027386226",
        "9728986075621437350131504894128984146939551938810073671231633620616345344412",
        "10579382052089733216628873394134968879891026686695240299956972154694558493896",
        "8171994519466002143995890536756742287314780571933910736618431096190430536601",
        "30420144259409274775063072923609924427757612539094840146996944760708902708570",
        "63962155989812703023698320394024694856871261481871757094333286947755599007133",
        "25280070391177856032024336895094721131222985610587247589336316615596140400436",
        "15305872319988027006162258914083163651002306183917888172691618513722838997098",
        "51545603291342006705870081001071419395633279951502747769141857387796043104608",
        "91109680756552587805002537489407348773333405839144382221272597323798859182191",
        "72175452855185658158184807496160149169667221240389196996344579971523681433202",
        "30361989157454953234766224747536334157139256334148153290771332849307087761025",
        "38169634499980959088614671703639492517637815232220682121652135514105493936992",
        "49591153263237620796156788742811547511792615129981565620486914545749079774827",
        "47403873018260745456113868791119169163627014766514972598212646481717066065016",
        "93989849689047144228924801010853106857960399638657695410345207191739048300111",
        "10590240512802509131776989274411792739339398409955259174829387591089799115255",
        "29183703335869638067547208413224742887766212046438654772943025958628178245227",
        "4131650227136944095885036960767735080970262672750406866066212532739784907379",
        "43395510588213653537697670365796375057855260611965666448183946252832290017444",
        "95246795133940226900907730059125298420936467652619708443128629427116119621152",
        "6012209003558496814495903476753006089125143165365334812097313083703216071080",
        "26183233284429251459198269925441295879550203824094631575778521083706115817955",
        "26058994700533582730528567480051558438548299522338811756875396252016497202713",
        "107240485663145290290374164860301805857261278222480421976433215167444496066511",
        "84412820763898503096477800002865877536719992495674955119188074297975154406587",
        "52386303852182662900790700046090769869460994629239741773176060026198900130384",
        "95746062835936512160025091603469309809932540674474329021370075533568318932379",
        "22711334660013961010382652754865456251782349529764119853461446587583972054666",
        "16959835233095757670013367728627149851239789174357906293937455553277911805495",
        "15116421110200928832147360650392633091242147433006813656250997138988179879750",
        "107878787525302837370688492081178689950008165750500003692400517211520334656293",
        "44210105558575948369921579518078229089923760124167628288943900602376706136436",
        "90305995748749060889452130219544332384396626628663475498252761213618628372367",
        "104941997925797907872686462815914481945432760720471803254797908465921520138024",
        "100036855232527386145662094141100441220151775745916101660987264242446845728894",
        "103285582836474146806606752170525767341430483568396209591447274936228630298052",
        "82197692939371228160449741709034077803239992888716859217989995857278406253737",
        "10040764964044995095453717286623030376397745892179877153575434454090155545240",
        "27304226040425863042893623786832369758179176309230053449707879364285977952630",
        "42627232144930751842910170221862679057276668485045156742021958050665662768084",
        "76972394926916659428228833084621905890924612368412796262119501852346293848159",
        "39796921406297542196667238133893946368231540421737718098283349901435707131075",
        "14745047092916651495052563068083093689676472592445845983334785004125684263162",
        "43421479365783318841667739359312715738029447177150400204380817518608837765863",
        "107871536756946365977710326147511195471121248998432910212631960353348700694610",
        "39505942243687894211614489736115535754716239859353578295470352855493707198619",
        "59676442091621150164811367352362126934419932715789994860508056194143441226580",
        "94470526851498636320865653968033227263836954414283116133326109455334870036212",
        "15044796858044094866329855531761112645684343559112419720568996573556805975600",
        "67157729293641241473980125231288476062565688273917759533572275886277269201651",
        "72911083146182058225942884942982388217243826839805061121973109250798137784134",
        "102973386186208530972563015865701244407271836208547629083437627219683649557477",
        "57485522356347377122696081086816661784954498123948319434326439317393351620564",
        "23112275556906805064863694321486306070917598599342299357379251070160695202292",
        "107618884362423342584703700349292347754139538760798319916678240538294838342400",
        "83961260400031958812820990908241261093246389047082613562825737834833753517337",
        "42726953951733266282750892844947149703751388034177248277671157488506520215317",
        "39379570934119946602507737250800178347029772561352879974941214627084076473292",
        "72203650529122342092280763801468513707870760755235719535090090101623606334441",
        "13389660788942842724553143053013919883368472759564135119390935439369513690496",
        "101745263541280877725997503552978999350831489463993178838531539940805924817361",
        "76849182334465191824607032600721023780793694103840553871800174717760598910761",
        "53896256317996838683363773836826653859512780625932638736752563553878867538095",
        "24688792501674999263943657175455335814404948006469220532686392550824320454904",
        "69132683906595821927803530074656979217668636557563597358799899743174233903941",
        "2861982085506615225917620192781928414994576134281371548401916333754363567986",
        "37311353286221616083824584705974993449107063556724405440534160586561042968316",
        "83718085796857523832195255218519255973031752296424117786202083986118546906913",
        "103633177691684814414226251117070754499104739002759424774194851613917008856616",
        "84968411062305024171594435878414659990735518025357685215223731503921265946461",
        "41865099330909055069143724769818364262362915440371474104937435863183989905059",
        "46156624920251322979270606388518884047396423747179340919303543598300663968593",
        "64416327466854458915398302811825971539792429791049027619471115285308743811583",
        "94942471312481523091911417289540395651121558150571128515230470225155209280585",
        "109682618775735319282534546194470743032129102295907200313471041846112653687024",
        "61531999191737540795124202104235799899980935519651613893518293245268304980543",
        "17797352534596268622733030076742840951214734697361029060619245779495726996632",
        "2323150752778983462106829021155031678603044899339819935981200101818542000989",
        "16998018904363448507967526489917057882529252665835717172712095240271574074587",
        "110634872413902251217040490777568744431854972018399530234679399294372694506842",
        "31639545145649753705216327198217551838008233610574104460826956396569310697060",
        "107845103764339268987018917144483935480716224058844669233389185480836000033760",
        "46240297572174662698030819651333060930818959915797061274854448535534474175039",
        "53065607123105696930220421963755520777674094852857308823370049733888025985616",
        "16931881300947470270453776207625163368485560075525342751440832370220475352149",
        "79254110800481916763656344422402393723573490114487681345184594841431920461089",
        "42268569642639492314994307446626647824927989776691987788682655102426770655233",
        "9749633319307409894058984489496091535125232227316143918000642155415596066903",
        "57606597628648270579042266322415267200058617178318601782866227410456726724976",
        "56082250485913115488341301630850455009935943641292622301678990296508134206571",
        "17957245764842844288802777667800779232762688847417238921175068882796163705248",
        "94356229516444419318132697346021621194464273500135725160277725602263001442644",
        "52536631226748676066386651084538409050048707922045928887930261833545619358914",
        "107794922118166328243581272159394479176678094739027519706768813902978100436849",
        "92984368734102511759118281503078145182557799453616537383408606074187034371208",
        "59652553897137603386525572460411404882917571255327541516871354737502335133690",
        "49012645345644326995052653072578673750516379655033952006121214224119466671764",
        "79025576845143484310735291619293962982804358365838961371044480386743856799994",
        "5437377540613244374799729812489584777222423091155743557287567155811057717409",
        "100687592213090267900708728796310211082532607828753010566886681655775031329660",
        "99074462968857696481475128596339544396152341206708424767062829343406495063192",
        "67476872698289965626550204192782761730653024363949045140720348870736942130242",
        "103307125141718054130755829916960708430672826104789971350239945481960770107890",
        "74087383014714668160537499936376991041273055222568604413015844459913259357334",
        "40924049099780965904051946083599822761993164889139026432053420731164022206736",
        "32594924940463736641240515015317856157169105212942308502676422036626316673214",
        "98990663138035055774586216545398054668349058134877723031747421828753359974443",
        "55821766022768786066770462759796825978667805772707620106340033118519147871694",
        "4001942224536365489828915551180230767516454384395893814399938353050969198154",
        "30136373426492646221252150708518703998248891683881870400906269276900707426865",
        "34943205764464817266133164313915763122699935186597909347522822673832250079664",
        "27737330737483170511275902246508559278973986181590368845166383812793468814968",
        "96292398813565494438359802278723334615526914389306923046282571355958508916558",
        "97147334956505986101750230325438660094766812949748276042292963837380833668274",
        "24754519562402723848413674701792328284127274989440581643644298347747941238812",
        "76111103490248669364580390783887028636436246028943665707064153006971943621186",
        "33764090322658516047637223655525551979364055499647855895233821795694749902854",
        "100536990630540359004783976190215234627391515555181073681294901127179838732969",
        "55991997435987096996680289872758998763908676069536901375395297778729059185671",
        "32860959903680178324832991459746631238726690317249285658471597044247794502256",
        "70074816806976994707467706079200635184034023598764203123459335544110485476930",
        "46213940675829172331116620705134022102338250410334045747023950259088879662946",
        "77000624259024986585504351395777746568094934279771127334532438603183524642061",
        "21719649576090832101273013788716623377603297433777804572370785470329817725170",
        "29209622978540575483991966565508890231057362045066230397327380085945876837821",
        "2445742484263083651472035320255578071935687960412507452207899496253120999364",
        "86846812580007547526361109808384103509272544750564766849178767957571523649544",
        "43025640639926253696325070988523609146060819319830735794100778654425057363895",
        "108957662689228031021948854644435971168708642184764962508575441689859324862868",
        "83891545396650121758556392255189778590486277180642660527000882403085396114823",
        "42527013475786190604202451803064937203698027000671529418992521798122995373551",
        "115180194520889678365425151865713593680657747284471744934804370945935167043862",
        "28979598171177052880917135045920701144584888536299261666846302083645491369348",
        "68351312608110279019109436395199010412431777911149851157132527077210966351650",
        "61759623963943995967580147094342313397376358019837276043205235302342147116585",
        "80714625408576660514217469096827255752431164791924432025682445176737446783085",
        "33048555646676368266608424610100449208381357250300222636992099726804869416731",
        "50682223610667325089810868083131721901859473966415125289975106060759036109476",
        "4271213571706787092297985431667190050727614825584809797590204884727103716461",
        "101314046722405990971733763321368296660561930294000591067108115987088407142646",
        "55565500177602146197728150332647093173137211885612327122425918553270191254877",
        "65556764608648687291293889343854786421750589271167654521933267288313526422497",
        "66877533773422945979143954094644173219583178339199697252673545117318799706373",
        "30511098623357801425494143655999121699575856091238269679669864984061501512835",
        "95900192636363991637086954986559552472749485926252879461208179855482821976623",
        "37879946127489462347049192209554168578320892231852882971030128420645686965013",
        "80479504274334215471057938992198620419540634144266821121799003865782336406529",
        "13326262422954139210095783388743602482455840337093117010479445267213907605425",
        "16047106134611124637925332265703907202779549268127518502853950466090054176776",
        "71499356105233640605079063493613576024353801558965221134519779175477723594865",
        "28438981751956157476540225984733791304599172905715743025543841239013139121102",
        "56066317647068426981453448715118237747130321302262827290362392918472904421147",
    ]
    .map(|y| F::from(BigUint::from_str(y).unwrap()))
    .to_vec();

    PoseidonConstants {
        round_keys,
        mds_matrix,
        num_full_rounds,
        num_partial_rounds,
    }
}
