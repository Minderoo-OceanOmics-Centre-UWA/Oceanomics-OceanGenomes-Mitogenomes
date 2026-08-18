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
