use ecfft::GoodCurve;

type Fp = halo2curves::secp256k1::Fp;

const CURVE_4_A: Fp = Fp::from_raw([
    1924362692430828527,
    180888387886949819,
    14444912836850558493,
    2716763698990320170,
]);

const CURVE_4_B_SQRT: Fp = Fp::from_raw([
    10596826214460559417,
    9041891995856355984,
    392200829566232436,
    5616829616257048236,
]);

const CURVE_4_GX: Fp = Fp::from_raw([
    3060553808241114122,
    4367422627483541323,
    1326591990371471461,
    1051615568340430255,
]);

const CURVE_4_GY: Fp = Fp::from_raw([
    1576479964359531032,
    10706990284747222844,
    2069836301523772900,
    11540652371418823164,
]);

const CURVE_4_CX: Fp = Fp::from_raw([
    1394469693679244729,
    3481743377114570646,
    685293755929734561,
    9752242693766949385,
]);

const CURVE_4_CY: Fp = Fp::from_raw([
    11112828892610998404,
    11816693849252775007,
    3142482327686601672,
    2138128838908646944,
]);

const CURVE_5_A: Fp = Fp::from_raw([
    18402892062958705657,
    10955586493449255806,
    274692491874833279,
    3521647190010012104,
]);

const CURVE_5_B_SQRT: Fp = Fp::from_raw([
    13277815540701934041,
    10000316683343802069,
    13748514902267845339,
    5043980866827216326,
]);

const CURVE_5_GX: Fp = Fp::from_raw([
    8681597433860724212,
    16010850653546434744,
    1655308633092053542,
    13482638234089226570,
]);

const CURVE_5_GY: Fp = Fp::from_raw([
    695535688134352662,
    12810977243071276429,
    6639318313449462421,
    9854099205183828948,
]);

const CURVE_5_CX: Fp = Fp::from_raw([
    7244846058583153822,
    15236867482366246868,
    7610066066648153412,
    8717324474930230203,
]);

const CURVE_5_CY: Fp = Fp::from_raw([
    15955524643521385563,
    14108119042026331605,
    8376852394104379031,
    5145942493709290957,
]);

const CURVE_6_A: Fp = Fp::from_raw([
    11754870036548954207,
    1758746815041297131,
    5040922207106606105,
    6156268686419792864,
]);

const CURVE_6_B_SQRT: Fp = Fp::from_raw([
    16551703456907310471,
    7307795367003411231,
    9107551177630293136,
    3643865576794489637,
]);

const CURVE_6_GX: Fp = Fp::from_raw([
    3057786712414561431,
    3924976238282577064,
    1535938406046208114,
    4471499328874959330,
]);

const CURVE_6_GY: Fp = Fp::from_raw([
    7420330572426678478,
    11093910355894798679,
    8046171174582240023,
    16159208434053522767,
]);

const CURVE_6_CX: Fp = Fp::from_raw([
    15232142892107882662,
    2997925312254061635,
    875684261157844424,
    8054980201271915862,
]);

const CURVE_6_CY: Fp = Fp::from_raw([
    5573271252396838460,
    7659129927758801858,
    11224891608690076565,
    8114225763096549468,
]);

const CURVE_7_A: Fp = Fp::from_raw([
    15267998901538414419,
    17985868627099147199,
    5570198032670981398,
    7365202425498739811,
]);

const CURVE_7_B_SQRT: Fp = Fp::from_raw([
    13238569970078865336,
    1859729155619525190,
    2289004025597154627,
    16424324367845100069,
]);

const CURVE_7_GX: Fp = Fp::from_raw([
    13067803014914932854,
    8460655374139991694,
    17522879348989963876,
    2592776320942502074,
]);

const CURVE_7_GY: Fp = Fp::from_raw([
    17581244616257969879,
    13563260062750024799,
    17836667944921387338,
    5158385585024810784,
]);

const CURVE_7_CX: Fp = Fp::from_raw([
    17627362889681060942,
    10449394617197091758,
    11211669951719111062,
    18402164978442722259,
]);

const CURVE_7_CY: Fp = Fp::from_raw([
    3635904622687808257,
    12660024001564793695,
    2997578841449112866,
    7489869964282615463,
]);

const CURVE_8_A: Fp = Fp::from_raw([
    16479441517948017563,
    12244661565122532810,
    16423402461885171455,
    15804938408404708752,
]);

const CURVE_8_B_SQRT: Fp = Fp::from_raw([
    4114471407724985276,
    6429895848762172356,
    9060307719139806083,
    1606308100763345976,
]);

const CURVE_8_GX: Fp = Fp::from_raw([
    18005553174453754936,
    7879246565041753863,
    15708703128473390087,
    12948592289805182905,
]);

const CURVE_8_GY: Fp = Fp::from_raw([
    2637815016833021192,
    5625620963822185667,
    15498097759340857613,
    2802364189360038003,
]);

const CURVE_8_CX: Fp = Fp::from_raw([
    12514982531648064548,
    7254771947927897203,
    6879061275311364813,
    4385541459413917142,
]);

const CURVE_8_CY: Fp = Fp::from_raw([
    13726278170638118925,
    10016993418218833106,
    13091102901943378213,
    8612533232618193985,
]);

const CURVE_9_A: Fp = Fp::from_raw([
    2821731813563793393,
    3977895281010832865,
    8603743292399951036,
    4645234720790204102,
]);

const CURVE_9_B_SQRT: Fp = Fp::from_raw([
    15890535715675950137,
    7339610358409226035,
    12609222160720627891,
    12499110658591842997,
]);

const CURVE_9_GX: Fp = Fp::from_raw([
    6103380741459351369,
    14746101125474414882,
    12417547802268400852,
    7335532149994146446,
]);

const CURVE_9_GY: Fp = Fp::from_raw([
    4181331351064768648,
    16489913493464340135,
    7051826832725726336,
    887431923330984487,
]);

const CURVE_9_CX: Fp = Fp::from_raw([
    9670988472649099633,
    15261760137634294840,
    2288914830631271678,
    6241984859397428357,
]);

const CURVE_9_CY: Fp = Fp::from_raw([
    3996701096097868069,
    16808707541849580191,
    2008740307070264540,
    10234541905633632584,
]);

const CURVE_10_A: Fp = Fp::from_raw([
    13443933661892288238,
    6366097774645666914,
    12539700524489124232,
    2960403700358460234,
]);

const CURVE_10_B_SQRT: Fp = Fp::from_raw([
    11179334656770650694,
    12204828351656968056,
    17469374953427230415,
    2698602761568343027,
]);

const CURVE_10_GX: Fp = Fp::from_raw([
    4752429915723436981,
    6658961595441054005,
    943316193080952835,
    10509103062531384873,
]);

const CURVE_10_GY: Fp = Fp::from_raw([
    7405820527339739030,
    1149755149636620515,
    12315441721581649311,
    9740641083146831387,
]);

const CURVE_10_CX: Fp = Fp::from_raw([
    12133915190353440166,
    12735241419273571667,
    984598181344074714,
    4945074058633718103,
]);

const CURVE_10_CY: Fp = Fp::from_raw([
    7787944055361603336,
    16188630343349344241,
    1798488611520969499,
    15905180573830923441,
]);

const CURVE_11_A: Fp = Fp::from_raw([
    44690967250983077,
    13024355091469571869,
    2426866618505792061,
    5439410159441159777,
]);

const CURVE_11_B_SQRT: Fp = Fp::from_raw([
    2482839174035592440,
    13977599229562359858,
    9165253311652858048,
    11796280965050311461,
]);

const CURVE_11_GX: Fp = Fp::from_raw([
    3785100838262116535,
    14366163517008314631,
    6520093107874784461,
    1432940500835404998,
]);

const CURVE_11_GY: Fp = Fp::from_raw([
    15446934078954168044,
    13724149936204307181,
    291296515805666972,
    17295416299404581082,
]);

const CURVE_11_CX: Fp = Fp::from_raw([
    3987179730389290606,
    12765099312359453542,
    14085665078244679772,
    1158541383839945849,
]);

const CURVE_11_CY: Fp = Fp::from_raw([
    2812404283588715887,
    10748530967036022352,
    15279323815639380689,
    7472866256744067949,
]);

pub fn secp256k1_good_curve(k: usize) -> (GoodCurve<Fp>, (Fp, Fp)) {
    if k == 4 {
        (
            GoodCurve::new(CURVE_4_A, CURVE_4_B_SQRT, CURVE_4_GX, CURVE_4_GY, k),
            (CURVE_4_CX, CURVE_4_CY),
        )
    } else if k == 5 {
        (
            GoodCurve::new(CURVE_5_A, CURVE_5_B_SQRT, CURVE_5_GX, CURVE_5_GY, k),
            (CURVE_5_CX, CURVE_5_CY),
        )
    } else if k == 6 {
        (
            GoodCurve::new(CURVE_6_A, CURVE_6_B_SQRT, CURVE_6_GX, CURVE_6_GY, k),
            (CURVE_6_CX, CURVE_6_CY),
        )
    } else if k == 7 {
        (
            GoodCurve::new(CURVE_7_A, CURVE_7_B_SQRT, CURVE_7_GX, CURVE_7_GY, k),
            (CURVE_7_CX, CURVE_7_CY),
        )
    } else if k == 8 {
        (
            GoodCurve::new(CURVE_8_A, CURVE_8_B_SQRT, CURVE_8_GX, CURVE_8_GY, k),
            (CURVE_8_CX, CURVE_8_CY),
        )
    } else if k == 9 {
        (
            GoodCurve::new(CURVE_9_A, CURVE_9_B_SQRT, CURVE_9_GX, CURVE_9_GY, k),
            (CURVE_9_CX, CURVE_9_CY),
        )
    } else if k == 10 {
        (
            GoodCurve::new(CURVE_10_A, CURVE_10_B_SQRT, CURVE_10_GX, CURVE_10_GY, k),
            (CURVE_10_CX, CURVE_10_CY),
        )
    } else if k == 11 {
        (
            GoodCurve::new(CURVE_11_A, CURVE_11_B_SQRT, CURVE_11_GX, CURVE_11_GY, k),
            (CURVE_11_CX, CURVE_11_CY),
        )
    } else {
        panic!("k must be between 4 and 11")
    }
}

#[cfg(test)]
mod tests {
    use ecfft::{find_coset_offset, GoodCurve};

    type F = halo2curves::secp256k1::Fp;

    fn to_limbs(x: F) -> [u64; 4] {
        let bytes = x.to_bytes();
        let mut limbs = [0u64; 4];

        for i in 0..4 {
            let mut limb_i = 0;
            for j in 0..8 {
                limb_i += (bytes[8 * i + j] as u64) << (8 * j);
            }
            limbs[i] = limb_i;
        }
        limbs
    }

    #[test]
    fn find_curves() {
        // We expect the tensor-IOP to use a square matrix for now,
        // so we only need to find curves with the square of the number
        // of evaluations
        for k in 4..12 {
            let curve = GoodCurve::<F>::find_k(k);
            let (coset_offset_x, coset_offset_y) =
                find_coset_offset(curve.a, curve.B_sqrt.square());
            println!(
                "const CURVE_{}_A: Fp = Fp::from_raw(
                {:?},
                );

                const CURVE_{}_B_SQRT: Fp = Fp::from_raw(
                    {:?},
                    );

                const CURVE_{}_GX: Fp = Fp::from_raw(
                    {:?},
                    );

                const CURVE_{}_GY: Fp = Fp::from_raw(
                {:?}
                );

                const CURVE_{}_CX: Fp = Fp::from_raw(
                    {:?},
                    );

                const CURVE_{}_CY: Fp = Fp::from_raw(
                {:?},
                );
                ",
                k,
                to_limbs(curve.a),
                k,
                to_limbs(curve.B_sqrt),
                k,
                to_limbs(curve.gx),
                k,
                to_limbs(curve.gy),
                k,
                to_limbs(coset_offset_x),
                k,
                to_limbs(coset_offset_y),
            );
        }
    }
}
