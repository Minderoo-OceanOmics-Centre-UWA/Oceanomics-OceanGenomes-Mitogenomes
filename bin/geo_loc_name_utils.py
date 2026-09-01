#!/usr/bin/env python3
"""
Shared geo_loc_name helpers.

Kept dependency-free so any script in bin/ can import it: Nextflow bind-mounts
the whole bin/ directory into the task container and puts it on PATH, so a
sibling import resolves via sys.path[0]. Mirrors species_name_utils.py.

geo_loc_name (the 'country' source modifier) is one of only two mandatory fields
in ENA's default sample checklist ERC000011, and its leading token must come from
the INSDC controlled vocabulary. The OceanOmics sample table is regenerated from a
spreadsheet, so values are corrected here on the way out rather than in the
database, where a fix would be undone by the next refresh.
"""


# INSDC geo_loc_name controlled vocabulary, taken from the
# geographic_location_country_andor_sea field of ENA checklist ERC000011:
#   https://www.ebi.ac.uk/ena/browser/api/xml/ERC000011
# Retrieved 2026-08-18. Re-check after an INSDC vocabulary update; the list is
# stable but not frozen.
INSDC_GEO_LOC_NAMES = frozenset({
    'Afghanistan', 'Albania', 'Algeria', 'American Samoa', 'Andorra', 'Angola', 'Anguilla',
    'Antarctica', 'Antigua and Barbuda', 'Arctic Ocean', 'Argentina', 'Armenia', 'Aruba',
    'Ashmore and Cartier Islands', 'Atlantic Ocean', 'Australia', 'Austria', 'Azerbaijan',
    'Bahamas', 'Bahrain', 'Baker Island', 'Baltic Sea', 'Bangladesh', 'Barbados',
    'Bassas da India', 'Belarus', 'Belgium', 'Belize', 'Benin', 'Bermuda', 'Bhutan',
    'Bolivia', 'Borneo', 'Bosnia and Herzegovina', 'Botswana', 'Bouvet Island', 'Brazil',
    'British Virgin Islands', 'Brunei', 'Bulgaria', 'Burkina Faso', 'Burundi', 'Cambodia',
    'Cameroon', 'Canada', 'Cape Verde', 'Cayman Islands', 'Central African Republic',
    'Chad', 'Chile', 'China', 'Christmas Island', 'Clipperton Island', 'Cocos Islands',
    'Colombia', 'Comoros', 'Cook Islands', 'Coral Sea Islands', 'Costa Rica',
    "Cote d'Ivoire", 'Croatia', 'Cuba', 'Curacao', 'Cyprus', 'Czech Republic', 'Czechia',
    'Democratic Republic of the Congo', 'Denmark', 'Djibouti', 'Dominica',
    'Dominican Republic', 'Ecuador', 'Egypt', 'El Salvador', 'Equatorial Guinea', 'Eritrea',
    'Estonia', 'Eswatini', 'Ethiopia', 'Europa Island', 'Falkland Islands (Islas Malvinas)',
    'Faroe Islands', 'Fiji', 'Finland', 'France', 'French Guiana', 'French Polynesia',
    'French Southern and Antarctic Lands', 'Gabon', 'Gambia', 'Gaza Strip', 'Georgia',
    'Germany', 'Ghana', 'Gibraltar', 'Glorioso Islands', 'Greece', 'Greenland', 'Grenada',
    'Guadeloupe', 'Guam', 'Guatemala', 'Guernsey', 'Guinea', 'Guinea-Bissau', 'Guyana',
    'Haiti', 'Heard Island and McDonald Islands', 'Honduras', 'Hong Kong', 'Howland Island',
    'Hungary', 'Iceland', 'India', 'Indian Ocean', 'Indonesia', 'Iran', 'Iraq', 'Ireland',
    'Isle of Man', 'Israel', 'Italy', 'Jamaica', 'Jan Mayen', 'Japan', 'Jarvis Island',
    'Jersey', 'Johnston Atoll', 'Jordan', 'Juan de Nova Island', 'Kazakhstan', 'Kenya',
    'Kerguelen Archipelago', 'Kingman Reef', 'Kiribati', 'Kosovo', 'Kuwait', 'Kyrgyzstan',
    'Laos', 'Latvia', 'Lebanon', 'Lesotho', 'Liberia', 'Libya', 'Liechtenstein',
    'Line Islands', 'Lithuania', 'Luxembourg', 'Macau', 'Madagascar', 'Malawi', 'Malaysia',
    'Maldives', 'Mali', 'Malta', 'Marshall Islands', 'Martinique', 'Mauritania',
    'Mauritius', 'Mayotte', 'Mediterranean Sea', 'Mexico',
    'Micronesia, Federated States of', 'Midway Islands', 'Moldova', 'Monaco', 'Mongolia',
    'Montenegro', 'Montserrat', 'Morocco', 'Mozambique', 'Myanmar', 'Namibia', 'Nauru',
    'Navassa Island', 'Nepal', 'Netherlands', 'New Caledonia', 'New Zealand', 'Nicaragua',
    'Niger', 'Nigeria', 'Niue', 'Norfolk Island', 'North Korea', 'North Macedonia',
    'North Sea', 'Northern Mariana Islands', 'Norway', 'Oman', 'Pacific Ocean', 'Pakistan',
    'Palau', 'Palmyra Atoll', 'Panama', 'Papua New Guinea', 'Paracel Islands', 'Paraguay',
    'Peru', 'Philippines', 'Pitcairn Islands', 'Poland', 'Portugal', 'Puerto Rico', 'Qatar',
    'Republic of the Congo', 'Reunion', 'Romania', 'Ross Sea', 'Russia', 'Rwanda',
    'Saint Barthelemy', 'Saint Helena', 'Saint Kitts and Nevis', 'Saint Lucia',
    'Saint Martin', 'Saint Pierre and Miquelon', 'Saint Vincent and the Grenadines',
    'Samoa', 'San Marino', 'Sao Tome and Principe', 'Saudi Arabia', 'Senegal', 'Serbia',
    'Seychelles', 'Sierra Leone', 'Singapore', 'Sint Maarten', 'Slovakia', 'Slovenia',
    'Solomon Islands', 'Somalia', 'South Africa',
    'South Georgia and the South Sandwich Islands', 'South Korea', 'South Sudan',
    'Southern Ocean', 'Spain', 'Spratly Islands', 'Sri Lanka', 'State of Palestine',
    'Sudan', 'Suriname', 'Svalbard', 'Sweden', 'Switzerland', 'Syria', 'Taiwan',
    'Tajikistan', 'Tanzania', 'Tasman Sea', 'Thailand', 'Timor-Leste', 'Togo', 'Tokelau',
    'Tonga', 'Trinidad and Tobago', 'Tromelin Island', 'Tunisia', 'Turkey', 'Turkmenistan',
    'Turks and Caicos Islands', 'Tuvalu', 'USA', 'Uganda', 'Ukraine',
    'United Arab Emirates', 'United Kingdom', 'Uruguay', 'Uzbekistan', 'Vanuatu',
    'Venezuela', 'Viet Nam', 'Virgin Islands', 'Wake Island', 'Wallis and Futuna',
    'West Bank', 'Western Sahara', 'Yemen', 'Zambia', 'Zimbabwe', 'missing',
    'missing: control sample', 'missing: data agreement established pre-2023',
    'missing: endangered species', 'missing: human-identifiable', 'missing: lab stock',
    'missing: sample group', 'missing: synthetic construct', 'missing: third party data',
    'not applicable', 'not collected', 'not provided', 'restricted access',
})


# Values the sample table holds that are a typo, a casing slip, or a long-form
# political name for a real INSDC country. Keyed lower-case; the replacement is
# written verbatim. Row counts are as of 2026-08-18.
GEO_LOC_ALIASES = {
    "republic of palau": "Palau",                       # 64 rows
    "japan": "Japan",                                   # 25 rows recorded as JAPAN
    "austalia": "Australia",                            # 15 rows, typo
    "kiribat": "Kiribati",                              # 3 rows, typo
    "falklands": "Falkland Islands (Islas Malvinas)",   # 1 row
    "kingdom of tonga": "Tonga",                        # 1 row
}


# Values that are not a country at all, so no single alias can be correct in
# general. Each one in the table today corresponds to a single survey, so it
# resolves on the country token plus the leading state token, and any OTHER use
# of the same country token is reported unmapped rather than guessed at. A None
# state means the country token alone already names one cruise.
# Oceans were read off the recorded coordinates, not assumed from the label.
SURVEY_GEO_LOC = {
    # 101 rows. lat -12.6..-14.0, lon 108.9..112.0, eastern Indian Ocean off WA.
    # The 16 rows with no coordinates carry the same 'Roo Rise' state.
    ("high seas", "roo rise"): "Indian Ocean",
    # 7 rows. lat 12.6..20.5, lon 151.1..156.5, western Pacific.
    ("international waters - chinese research vessel", None): "Pacific Ocean",
    # 1 row at -61.140, -58.019, north of the Antarctic Peninsula. Marine
    # specimen; 'Antarctica' would be the value for a terrestrial collection.
    ("south shetland", None): "Southern Ocean",
}


# INSDC missing-value term for a mandatory field with no recorded value. Omitting
# geo_loc_name is not an option: it is one of only two mandatory fields in
# ERC000011, so a blank one fails sample registration outright. "not provided"
# states the absence explicitly, which is what the vocabulary provides it for.
MISSING_GEO_LOC = "not provided"


# Country recovered from an NCBI BioSample this project has already registered,
# for og_ids whose sample row carries no country and no locality. These records
# are INSDC metadata we submitted ourselves, so they are authoritative rather
# than inferred, and using them here beats writing a missing-value term.
#
# Built by fetching every ncbi_biosample_id belonging to a mitogenome sample with
# no country and no location, on 2026-08-18:
#   efetch.fcgi?db=biosample&id=<accessions>&rettype=xml
# Each record's geo_loc_name was kept only when its `isolate` attribute matched
# the og_id independently, so a mislinked accession cannot inject a wrong
# country; all 132 matched. No coordinates were recoverable (lat_lon is
# 'not provided' on every record). Re-derive with the same query if the set of
# registered BioSamples grows.
BIOSAMPLE_GEO_LOC = {
    "OG260": "Australia", "OG261": "Australia", "OG262": "Australia", "OG266": "Australia",
    "OG267": "Australia", "OG268": "Australia", "OG271": "Australia", "OG272": "Australia",
    "OG273": "Australia", "OG274": "Australia", "OG275": "Australia", "OG276": "Australia",
    "OG277": "Australia", "OG278": "Australia", "OG280": "Australia", "OG281": "Australia",
    "OG286": "Australia", "OG290": "Australia", "OG291": "Australia", "OG292": "Australia",
    "OG293": "Australia", "OG299": "Australia", "OG301": "Australia", "OG302": "Australia",
    "OG305": "Australia", "OG307": "Australia", "OG315": "Australia", "OG316": "Australia",
    "OG320": "Australia", "OG321": "Australia", "OG329": "Australia", "OG331": "Australia",
    "OG374": "Australia", "OG376": "Australia", "OG377": "Australia", "OG378": "Australia",
    "OG380": "Australia", "OG381": "Australia", "OG382": "Australia", "OG384": "Australia",
    "OG385": "Australia", "OG386": "Australia", "OG387": "Australia", "OG388": "Australia",
    "OG389": "Australia", "OG390": "Australia", "OG391": "Australia", "OG401": "Australia",
    "OG402": "Australia", "OG403": "Australia", "OG405": "Australia", "OG410": "Australia",
    "OG411": "Australia", "OG413": "Australia", "OG416": "Australia", "OG418": "Australia",
    "OG419": "Australia", "OG446": "Australia", "OG450": "Australia", "OG454": "Australia",
    "OG456": "Australia", "OG459": "Australia", "OG461": "Australia", "OG463": "Australia",
    "OG465": "Australia", "OG469": "Australia", "OG477": "Australia", "OG478": "Australia",
    "OG479": "Australia", "OG480": "Australia", "OG482": "Australia", "OG484": "Australia",
    "OG492": "Australia", "OG493": "Australia", "OG494": "Australia", "OG495": "Australia",
    "OG497": "Australia", "OG498": "Australia", "OG501": "Australia", "OG504": "Australia",
    "OG505": "Australia", "OG509": "Australia", "OG517": "Australia", "OG518": "Australia",
    "OG519": "Australia", "OG520": "Australia", "OG521": "Australia", "OG522": "Australia",
    "OG523": "Australia", "OG531": "Australia", "OG538": "Australia", "OG544": "Australia",
    "OG545": "Australia", "OG546": "Australia", "OG547": "Australia", "OG549": "Australia",
    "OG550": "Australia", "OG552": "Australia", "OG553": "Australia", "OG554": "Australia",
    "OG556": "Australia", "OG557": "Australia", "OG558": "Australia", "OG559": "Australia",
    "OG560": "Australia", "OG563": "Australia", "OG565": "Australia", "OG574": "Australia",
    "OG575": "Australia", "OG577": "Australia", "OG578": "Australia", "OG581": "Australia",
    "OG585": "Australia", "OG586": "Australia", "OG587": "Australia", "OG588": "Australia",
    "OG589": "Australia", "OG590": "Australia", "OG591": "Australia", "OG594": "Australia",
    "OG595": "Australia", "OG597": "Australia", "OG599": "Australia", "OG600": "Australia",
    "OG601": "Australia", "OG603": "Australia", "OG604": "Australia", "OG612": "Australia",
    "OG614": "Australia", "OG615": "Australia", "OG621": "Australia", "OG624": "Australia",
}


# Hemispheres implied by an INSDC geo_loc_name, as (latitude, longitude).
#
# The sample table records latitudes as unsigned magnitudes -- a Ningaloo sample
# at 22.03 S is stored as "22.03 113.891" -- so the hemisphere has to come from
# somewhere else. It comes from here. table2asn checks the coordinate against the
# country and raises SEQ_DESCR.LatLonValue when they disagree, which is exactly
# the check being anticipated.
#
# A None means the territory genuinely falls on both sides of the equator or of
# the prime/anti meridian, so no single answer is correct and the caller must
# omit the coordinate rather than pick one. A country absent from this table is
# treated the same way. That is deliberately the safe direction: an omitted
# lat_lon is recoverable, a confidently wrong one submitted to a public archive
# is not. Add entries as new collection localities appear.
#
# Bounds were taken from each territory's full extent including outlying islands,
# which is why several countries that read as unambiguous are None: New Zealand's
# Kermadec and Chatham Islands sit west of 180, Fiji and Kiribati and Russia and
# the Aleutians straddle the antimeridian, and the United Kingdom, France, Spain,
# Algeria, Mali, Burkina Faso and Ghana straddle the prime meridian.
COUNTRY_HEMISPHERE = {
    # Australia and the eastern Indian Ocean
    'Australia': ('S', 'E'), 'Ashmore and Cartier Islands': ('S', 'E'),
    'Christmas Island': ('S', 'E'), 'Cocos Islands': ('S', 'E'),
    'Coral Sea Islands': ('S', 'E'), 'Norfolk Island': ('S', 'E'),
    'Timor-Leste': ('S', 'E'), 'Indonesia': (None, 'E'), 'Borneo': (None, 'E'),
    'Malaysia': ('N', 'E'), 'Brunei': ('N', 'E'), 'Singapore': ('N', 'E'),
    'Papua New Guinea': ('S', 'E'), 'Solomon Islands': ('S', 'E'),
    'Vanuatu': ('S', 'E'), 'New Caledonia': ('S', 'E'),
    'New Zealand': ('S', None), 'Fiji': ('S', None),

    # Pacific
    'Nauru': ('S', 'E'), 'Tuvalu': ('S', 'E'), 'Marshall Islands': ('N', 'E'),
    'Micronesia, Federated States of': ('N', 'E'), 'Palau': ('N', 'E'),
    'Guam': ('N', 'E'), 'Northern Mariana Islands': ('N', 'E'),
    'Wake Island': ('N', 'E'), 'Kiribati': (None, None),
    'Tonga': ('S', 'W'), 'Samoa': ('S', 'W'), 'American Samoa': ('S', 'W'),
    'Niue': ('S', 'W'), 'Tokelau': ('S', 'W'), 'Cook Islands': ('S', 'W'),
    'Wallis and Futuna': ('S', 'W'), 'French Polynesia': ('S', 'W'),
    'Pitcairn Islands': ('S', 'W'), 'Jarvis Island': ('S', 'W'),
    'Midway Islands': ('N', 'W'), 'Johnston Atoll': ('N', 'W'),
    'Palmyra Atoll': ('N', 'W'), 'Kingman Reef': ('N', 'W'),
    'Howland Island': ('N', 'W'), 'Baker Island': ('N', 'W'),
    'Clipperton Island': ('N', 'W'), 'Line Islands': (None, 'W'),
    'Spratly Islands': ('N', 'E'), 'Paracel Islands': ('N', 'E'),

    # Asia
    'Japan': ('N', 'E'), 'China': ('N', 'E'), 'Taiwan': ('N', 'E'),
    'Hong Kong': ('N', 'E'), 'Macau': ('N', 'E'), 'South Korea': ('N', 'E'),
    'North Korea': ('N', 'E'), 'Mongolia': ('N', 'E'), 'Philippines': ('N', 'E'),
    'Viet Nam': ('N', 'E'), 'Laos': ('N', 'E'), 'Cambodia': ('N', 'E'),
    'Thailand': ('N', 'E'), 'Myanmar': ('N', 'E'), 'Bangladesh': ('N', 'E'),
    'Bhutan': ('N', 'E'), 'Nepal': ('N', 'E'), 'India': ('N', 'E'),
    'Sri Lanka': ('N', 'E'), 'Maldives': (None, 'E'), 'Pakistan': ('N', 'E'),
    'Afghanistan': ('N', 'E'), 'Kazakhstan': ('N', 'E'), 'Uzbekistan': ('N', 'E'),
    'Turkmenistan': ('N', 'E'), 'Tajikistan': ('N', 'E'), 'Kyrgyzstan': ('N', 'E'),
    'Russia': ('N', None),

    # Middle East
    'Iran': ('N', 'E'), 'Iraq': ('N', 'E'), 'Kuwait': ('N', 'E'),
    'Saudi Arabia': ('N', 'E'), 'Bahrain': ('N', 'E'), 'Qatar': ('N', 'E'),
    'United Arab Emirates': ('N', 'E'), 'Oman': ('N', 'E'), 'Yemen': ('N', 'E'),
    'Jordan': ('N', 'E'), 'Israel': ('N', 'E'), 'Lebanon': ('N', 'E'),
    'Syria': ('N', 'E'), 'Turkey': ('N', 'E'), 'Cyprus': ('N', 'E'),
    'Gaza Strip': ('N', 'E'), 'West Bank': ('N', 'E'),
    'State of Palestine': ('N', 'E'), 'Georgia': ('N', 'E'),
    'Armenia': ('N', 'E'), 'Azerbaijan': ('N', 'E'),

    # Europe
    'Iceland': ('N', 'W'), 'Ireland': ('N', 'W'), 'United Kingdom': ('N', None),
    'Isle of Man': ('N', 'W'), 'Jersey': ('N', 'W'), 'Guernsey': ('N', 'W'),
    'Faroe Islands': ('N', 'W'), 'Greenland': ('N', 'W'), 'Jan Mayen': ('N', 'W'),
    'Svalbard': ('N', 'E'), 'Norway': ('N', 'E'), 'Sweden': ('N', 'E'),
    'Finland': ('N', 'E'), 'Denmark': ('N', 'E'), 'Estonia': ('N', 'E'),
    'Latvia': ('N', 'E'), 'Lithuania': ('N', 'E'), 'Belarus': ('N', 'E'),
    'Ukraine': ('N', 'E'), 'Moldova': ('N', 'E'), 'Poland': ('N', 'E'),
    'Germany': ('N', 'E'), 'Netherlands': ('N', 'E'), 'Belgium': ('N', 'E'),
    'Luxembourg': ('N', 'E'), 'France': ('N', None), 'Monaco': ('N', 'E'),
    'Andorra': ('N', 'E'), 'Spain': ('N', None), 'Gibraltar': ('N', 'W'),
    'Portugal': ('N', 'W'), 'Italy': ('N', 'E'), 'San Marino': ('N', 'E'),
    'Malta': ('N', 'E'), 'Switzerland': ('N', 'E'), 'Liechtenstein': ('N', 'E'),
    'Austria': ('N', 'E'), 'Czech Republic': ('N', 'E'), 'Czechia': ('N', 'E'),
    'Slovakia': ('N', 'E'), 'Hungary': ('N', 'E'), 'Slovenia': ('N', 'E'),
    'Croatia': ('N', 'E'), 'Bosnia and Herzegovina': ('N', 'E'),
    'Serbia': ('N', 'E'), 'Kosovo': ('N', 'E'), 'Montenegro': ('N', 'E'),
    'North Macedonia': ('N', 'E'), 'Albania': ('N', 'E'), 'Greece': ('N', 'E'),
    'Bulgaria': ('N', 'E'), 'Romania': ('N', 'E'),

    # Africa
    'Morocco': ('N', 'W'), 'Western Sahara': ('N', 'W'), 'Algeria': ('N', None),
    'Tunisia': ('N', 'E'), 'Libya': ('N', 'E'), 'Egypt': ('N', 'E'),
    'Sudan': ('N', 'E'), 'South Sudan': ('N', 'E'), 'Eritrea': ('N', 'E'),
    'Djibouti': ('N', 'E'), 'Ethiopia': ('N', 'E'), 'Somalia': (None, 'E'),
    'Kenya': (None, 'E'), 'Uganda': (None, 'E'), 'Rwanda': ('S', 'E'),
    'Burundi': ('S', 'E'), 'Tanzania': ('S', 'E'), 'Malawi': ('S', 'E'),
    'Zambia': ('S', 'E'), 'Zimbabwe': ('S', 'E'), 'Mozambique': ('S', 'E'),
    'Botswana': ('S', 'E'), 'Namibia': ('S', 'E'), 'South Africa': ('S', 'E'),
    'Lesotho': ('S', 'E'), 'Eswatini': ('S', 'E'), 'Angola': ('S', 'E'),
    'Democratic Republic of the Congo': (None, 'E'),
    'Republic of the Congo': (None, 'E'), 'Gabon': (None, 'E'),
    'Equatorial Guinea': (None, 'E'), 'Sao Tome and Principe': (None, 'E'),
    'Cameroon': ('N', 'E'), 'Central African Republic': ('N', 'E'),
    'Chad': ('N', 'E'), 'Niger': ('N', 'E'), 'Nigeria': ('N', 'E'),
    'Benin': ('N', 'E'), 'Togo': ('N', 'E'), 'Ghana': ('N', None),
    'Burkina Faso': ('N', None), 'Mali': ('N', None),
    "Cote d'Ivoire": ('N', 'W'), 'Liberia': ('N', 'W'),
    'Sierra Leone': ('N', 'W'), 'Guinea': ('N', 'W'),
    'Guinea-Bissau': ('N', 'W'), 'Senegal': ('N', 'W'), 'Gambia': ('N', 'W'),
    'Mauritania': ('N', 'W'), 'Cape Verde': ('N', 'W'),
    'Saint Helena': ('S', 'W'), 'Madagascar': ('S', 'E'), 'Comoros': ('S', 'E'),
    'Mayotte': ('S', 'E'), 'Mauritius': ('S', 'E'), 'Reunion': ('S', 'E'),
    'Seychelles': ('S', 'E'), 'Europa Island': ('S', 'E'),
    'Bassas da India': ('S', 'E'), 'Juan de Nova Island': ('S', 'E'),
    'Tromelin Island': ('S', 'E'), 'Glorioso Islands': ('S', 'E'),

    # Americas
    'Canada': ('N', 'W'), 'USA': ('N', None), 'Mexico': ('N', 'W'),
    'Guatemala': ('N', 'W'), 'Belize': ('N', 'W'), 'Honduras': ('N', 'W'),
    'El Salvador': ('N', 'W'), 'Nicaragua': ('N', 'W'), 'Costa Rica': ('N', 'W'),
    'Panama': ('N', 'W'), 'Cuba': ('N', 'W'), 'Jamaica': ('N', 'W'),
    'Haiti': ('N', 'W'), 'Dominican Republic': ('N', 'W'),
    'Puerto Rico': ('N', 'W'), 'Bahamas': ('N', 'W'), 'Bermuda': ('N', 'W'),
    'Cayman Islands': ('N', 'W'), 'Turks and Caicos Islands': ('N', 'W'),
    'British Virgin Islands': ('N', 'W'), 'Virgin Islands': ('N', 'W'),
    'Anguilla': ('N', 'W'), 'Antigua and Barbuda': ('N', 'W'),
    'Saint Kitts and Nevis': ('N', 'W'), 'Montserrat': ('N', 'W'),
    'Guadeloupe': ('N', 'W'), 'Dominica': ('N', 'W'), 'Martinique': ('N', 'W'),
    'Saint Lucia': ('N', 'W'), 'Saint Vincent and the Grenadines': ('N', 'W'),
    'Grenada': ('N', 'W'), 'Barbados': ('N', 'W'),
    'Trinidad and Tobago': ('N', 'W'), 'Aruba': ('N', 'W'), 'Curacao': ('N', 'W'),
    'Sint Maarten': ('N', 'W'), 'Saint Martin': ('N', 'W'),
    'Saint Barthelemy': ('N', 'W'), 'Saint Pierre and Miquelon': ('N', 'W'),
    'Navassa Island': ('N', 'W'), 'Colombia': (None, 'W'),
    'Venezuela': ('N', 'W'), 'Guyana': ('N', 'W'), 'Suriname': ('N', 'W'),
    'French Guiana': ('N', 'W'), 'Ecuador': (None, 'W'), 'Brazil': (None, 'W'),
    'Peru': ('S', 'W'), 'Bolivia': ('S', 'W'), 'Paraguay': ('S', 'W'),
    'Uruguay': ('S', 'W'), 'Argentina': ('S', 'W'), 'Chile': ('S', 'W'),
    'Falkland Islands (Islas Malvinas)': ('S', 'W'),
    'South Georgia and the South Sandwich Islands': ('S', 'W'),

    # Polar and marine
    'Antarctica': ('S', None), 'Bouvet Island': ('S', 'E'),
    'Heard Island and McDonald Islands': ('S', 'E'),
    'Kerguelen Archipelago': ('S', 'E'),
    'French Southern and Antarctic Lands': ('S', 'E'),
    'Arctic Ocean': ('N', None), 'Southern Ocean': ('S', None),
    'Ross Sea': ('S', None), 'Atlantic Ocean': (None, None),
    'Pacific Ocean': (None, None), 'Indian Ocean': (None, None),
    'North Sea': ('N', 'E'), 'Baltic Sea': ('N', 'E'),
    'Mediterranean Sea': ('N', None), 'Tasman Sea': ('S', 'E'),
}


def hemispheres_for_geo_loc_name(value):
    """Return the (latitude, longitude) hemispheres a geo_loc_name implies.

    `value` is a resolved geo_loc_name, "<country>[: <locality>]"; only the
    leading country token is consulted. Either element is None when the country
    straddles that axis, or is not in COUNTRY_HEMISPHERE at all -- including
    every INSDC missing-value term, which names no place. A caller with a
    coordinate that states no hemisphere of its own must omit it in that case
    rather than guess.
    """
    country = str(value or "").partition(":")[0].strip()
    return COUNTRY_HEMISPHERE.get(country, (None, None))


def _rejoin(country, remainder):
    """Reattach the locality suffix, if any, to a replaced country token."""
    return f"{country}: {remainder}" if remainder else country


def resolve_geo_loc_name(value, og_id=None):
    """Map a sample-table geo_loc_name onto the INSDC controlled vocabulary.

    `value` is the composite the metadata query builds, "<country>[: <state>][,
    <locality>]" (see the CASE expression in build_source_modifiers.py). Only the
    leading country token is controlled, so any locality suffix is carried
    through untouched. `og_id` is optional and used only to recover a country
    from an already-registered BioSample when the row records none.

    Returns (value, status) where status is one of:
        'missing'  -- nothing recorded anywhere; resolves to the INSDC
                      'not provided' term rather than an omitted (and therefore
                      invalid) field
        'biosample'-- nothing on the sample row, but an NCBI BioSample we already
                      registered for this og_id states the country
        'ok'       -- already a controlled value, returned unchanged
        'aliased'  -- rewritten onto a controlled value
        'derived'  -- no country column, but the locality named one
        'unmapped' -- not controlled and no alias, returned UNCHANGED so that
                      table2asn raises SEQ_DESCR.BadGeoLocNameCode and the
                      validation gate quarantines that one assembly. Substituting
                      a placeholder here would ship a sample with its locality
                      quietly replaced, which is worse than a visible failure.
    """
    text = str(value or "").strip()
    if not text or text.lower() in ("unknown", "original locality unknown"):
        # Nothing recorded on the sample row, but this project may already have
        # registered a BioSample for it that states the country. Consulted only
        # here, so a recorded country always wins over the recovered one.
        recovered = BIOSAMPLE_GEO_LOC.get(str(og_id or "").strip())
        if recovered:
            return recovered, "biosample"
        return MISSING_GEO_LOC, "missing"

    country, _sep, remainder = text.partition(":")
    country, remainder = country.strip(), remainder.strip()
    key = country.lower()

    for (survey_country, survey_state), resolved in SURVEY_GEO_LOC.items():
        if key == survey_country and (
            survey_state is None or remainder.lower().startswith(survey_state)
        ):
            return _rejoin(resolved, remainder), "aliased"

    alias = GEO_LOC_ALIASES.get(key)
    if alias:
        # A case-only alias (JAPAN -> Japan) still counts as a rewrite; an exact
        # match means the value was already correct.
        return _rejoin(alias, remainder), "ok" if alias == country else "aliased"

    if country in INSDC_GEO_LOC_NAMES:
        return text, "ok"

    # No country column, but the locality itself may name one: rows with a NULL
    # country reach here as their bare location text ("Israel, Elat, Gulf of
    # Aquaba"). Promote the leading comma-separated token when it is a controlled
    # value, keeping the rest as the locality. Only ever consulted for a value
    # that has already failed every lookup above, so it cannot override a
    # recorded country.
    if not remainder:
        head, _comma, tail = text.partition(",")
        head, tail = head.strip(), tail.strip()
        if head in INSDC_GEO_LOC_NAMES:
            return _rejoin(head, tail), "derived"

    return text, "unmapped"


def unmapped_warning(seqid, value):
    """Operator-facing message naming the exact line to add for a new value."""
    country = str(value or "").partition(":")[0].strip()
    return (
        f"⚠️  {seqid}: geo_loc_name {country!r} is not an INSDC controlled "
        f"value and has no alias, so this assembly will be quarantined.\n"
        f"    Add to bin/geo_loc_name_utils.py GEO_LOC_ALIASES:\n"
        f'        "{country.lower()}": "<correct INSDC value>",'
    )
