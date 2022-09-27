import json

problem = {
    "num_vars": 16,
    "names": [
        "bb",
        "satdk",
        "satpsi",
        "slop",
        "smcmax",
        "wltsmc",
        "D",
        "coeff_secondary",
        "exponent_secondary",
        "max_gw_storage",
        "Cgw",
        "expon",
        "K_nash",
        "refkdt",
        "trigger_z_m",
        "fc_atm_press_fraction"
    ],
    "bounds": [
        [
            10,
            10
        ],
        [
            0.00013888888,
            0.00013888888
        ],
        [
            0.657124897,
            0.657124897
        ],
        [
            1,
            1
        ],
        [
            0.418156524,
            0.418156524
        ],
        [
            0.310002781,
            0.310002781
        ],
        [
            0.87,
            0.87
        ],
        [
            0.5,
            0.5
        ],
        [
            1,
            1
        ],
        [
            0.5,
            0.5
        ],
        [
            5,
            5
        ],
        [
            0.0002,
            0.0002
        ],
        [
            0.1,
            0.1
        ],
        [
            0.496733933,
            0.496733933
        ],
        [
            0.109099529,
            0.109099529
        ],
        [
            0.33,
            0.33
        ]
    ]
}

fn = r'G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\2_data_input\config_json_SALibtemp.json'
with open(fn, 'w') as outfile:
    json.dump(problem, outfile, indent=4)
