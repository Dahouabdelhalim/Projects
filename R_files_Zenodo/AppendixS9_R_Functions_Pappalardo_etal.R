
# ------------------------------------------------------------------------
# -------Comparing methods for mapping global parasite diversity ---------
# ------------------------------------------------------------------------
#
# Most functions were written by Paula Pappalardo, the functions written by
# Ignacio Morales-Castilla are highlighted.


# ------Functions to filter GMPD dataset and fix coordinates problems
# -------------------------------------------------------------------
fixCoordinates <- function(mydf){
  # Fix geographic coordinates reported in GMPDv2 Stephens et al. (2017)
  # (a few random problems fixed, most of them points falling in water)
  #
  # Args:
  #   mydf: dataframe with columns "Longitude","Latitude", "Citation", "LocationName"
  # Returns:
  #   same dataframe with issues fixed
  # Identify lines with problems and fixed them
  edited.df <- mydf %>% dplyr::mutate(Latitude= ifelse(LocationName == "Okongena, Lubutu, Belgian Congo", -3.458261, Latitude),
                                      Longitude= ifelse(LocationName == "Okongena, Lubutu, Belgian Congo", 23.05678, Longitude)) %>% 
    dplyr::mutate(Latitude= ifelse(LocationName == "Rondo", 3.385755, Latitude),
                  Longitude= ifelse(LocationName == "Rondo", -76.53993, Longitude),
                  LocationName= ifelse(LocationName == "Rondo", "Rondo, Colombia", LocationName)) %>%
    dplyr::mutate(Latitude= ifelse(LocationName == "Tacheng", 24.47528, Latitude),
                  Longitude= ifelse(LocationName == "Tacheng", 101.3431, Longitude)) %>% 
    dplyr::mutate(Latitude= ifelse(LocationName == "Barbascal", 3.68334, Latitude),
                  Longitude= ifelse(LocationName == "Barbascal", -73.43331, Longitude),
                  LocationName= ifelse(LocationName == "Barbascal", "Barbascal, Colombia", LocationName)) %>% 
    dplyr::mutate(Latitude= ifelse(LocationName == "Pemba", -5.031935, Latitude),
                  Longitude= ifelse(LocationName == "Pemba", 39.77556, Longitude),
                  LocationName= ifelse(LocationName == "Pemba", "Pemba, Zanzibar", LocationName)) %>% 
    dplyr::mutate(Latitude= ifelse(LocationName == "west of Great Slave Lake in the Northwest territories", 61.441, Latitude),
                  Longitude= ifelse(LocationName == "west of Great Slave Lake in the Northwest territories", -116.823, Longitude),
                  LocationName= ifelse(LocationName == "west of Great Slave Lake in the Northwest territories", "west of Great Slave Lake, Northwest territories, Canada", LocationName)) %>% 
    dplyr::mutate(Latitude= ifelse(LocationName == "Iwate Prefecture, Japan", 39.704, Latitude),
                  Longitude= ifelse(LocationName == "Iwate Prefecture, Japan", 141.153, Longitude)) %>% 
    dplyr::mutate(Longitude= ifelse(LocationName == "Northern Lower Peninsula, Michigan (counties Montmorency and Otsego)", -84.14, Longitude)) %>% # it needed to be negative for West
    dplyr::mutate(Latitude= ifelse(LocationName == "North and Central Florida", 29.522, Latitude),
                  Longitude= ifelse(LocationName == "North and Central Florida", -82.041, Longitude)) %>% 
    dplyr::mutate(Latitude= ifelse(LocationName == "Northwestern Boehmia Sumny Dul, Ore Mountains, Czech Republic (Krusne hory in czech)", 50.58, Latitude),
                  Longitude= ifelse(LocationName == "Northwestern Boehmia Sumny Dul, Ore Mountains, Czech Republic (Krusne hory in czech)",  13, Longitude)) %>% 
    dplyr::mutate(Latitude= ifelse(LocationName == "Yokohama, JAPAN", 35.444, Latitude),
                  Longitude= ifelse(LocationName == "Yokohama, JAPAN",  139.638, Longitude)) %>%
    # Locations that were not georeferenced in GMPD (have no binomial name..)
    dplyr::mutate(Latitude= ifelse(LocationName == "Aubrey Valley, Arizona", 35.533333, Latitude),
                  Longitude= ifelse(LocationName == "Aubrey Valley, Arizona",  -113.166667, Longitude)) %>% # from original paper
    dplyr::mutate(Latitude= ifelse(LocationName == "Buffalo Gap National Grasslands and Badlands National Park, South Dakota", 43.25, Latitude),
                  Longitude= ifelse(LocationName == "Buffalo Gap National Grasslands and Badlands National Park, South Dakota",   -109.1, Longitude)) %>% # from original paper
    dplyr::mutate(Latitude= ifelse(LocationName == "Janos, MEXICO", 31.066667, Latitude),
                  Longitude= ifelse(LocationName == "Janos, MEXICO",   -107.85, Longitude)) %>% # from original paper
    dplyr::mutate(Latitude= ifelse(LocationName == "Coyote Basin, Utah and Wolf Creek, Colorado", 40.25, Latitude),
                  Longitude= ifelse(LocationName == "Coyote Basin, Utah and Wolf Creek, Colorado",   -109.1, Longitude)) %>% # from original paper
    dplyr::mutate(Latitude= ifelse(LocationName == "Shirley Basin, Wyoming", 42.116667, Latitude),
                  Longitude= ifelse(LocationName == "Shirley Basin, Wyoming",   -106.05, Longitude)) %>% # from original paper
    dplyr::mutate(Latitude= ifelse(LocationName == "Meeteetse, Wyoming", 44.15718, Latitude),
                  Longitude= ifelse(LocationName == "Meeteetse, Wyoming",   -108.8715, Longitude)) %>% # from geocode
    #
    # Locations that were geocoded (using the centroid) that can be improved
    # --------Njiokou et al 2004, Cameroon, paper says southern Cameroon and list 4 localities with coordinates
    # but doesn't specified localities by host and parasite; most field surveys in Bipindi, coords from paper
    dplyr::mutate(Latitude= ifelse(Citation == "Njiokou et al 2004", 3.033333, Latitude),
                  Longitude= ifelse(Citation == "Njiokou et al 2004",   10.366667, Longitude)) %>%
    # ------Hugot 1985: paper has more details for some locations than GMPD
    # l'Arataye is the locality of French Guiana for Trypanoxyuris trypanuris, coords from Wikipedia for the source
    dplyr::mutate(Latitude= ifelse(Citation == "Hugot 1985" & ParasiteCorrectedName == "Trypanoxyuris trypanuris", 3.741944, Latitude),
                  Longitude= ifelse(Citation == "Hugot 1985" & ParasiteCorrectedName == "Trypanoxyuris trypanuris", -53.167222, Longitude)) %>%
    # Carimagua, Colombia is the locality for Trypanoxyuris microon (is on Puerto Gaitan), geocoded in R
    dplyr::mutate(Latitude= ifelse(Citation == "Hugot 1985" & ParasiteCorrectedName == "Trypanoxyuris microon", 4.612817, Latitude),
                  Longitude= ifelse(Citation == "Hugot 1985" & ParasiteCorrectedName == "Trypanoxyuris microon", -74.14528, Longitude)) %>%
    # Updated geocode for Venezuela
    dplyr::mutate(Latitude= ifelse(LocationName == "Venezuela", 6.42375, Latitude),
                  Longitude= ifelse(LocationName == "Venezuela", -66.58973, Longitude)) %>%
    # Updated geocode for Brazil
    dplyr::mutate(Latitude= ifelse(LocationName == "Brazil", -14.235, Latitude),
                  Longitude= ifelse(LocationName == "Brazil", -51.92528, Longitude)) %>%
    # Updated geocode for Colombia (paper had the French spelling)
    dplyr::mutate(Latitude= ifelse(LocationName == "Columbia",4.570868, Latitude),
                  Longitude= ifelse(LocationName == "Columbia", -74.29733, Longitude)) %>%
    # Updated geocode for Guyana
    dplyr::mutate(Latitude= ifelse(LocationName == "Guyana", 3.933889, Latitude),
                  Longitude= ifelse(LocationName == "Guyana",-53.12578, Longitude)) %>%
    # Updated geocode for Roraima
    dplyr::mutate(Latitude= ifelse(LocationName == "Roraima", 2.737597, Latitude),
                  Longitude= ifelse(LocationName == "Roraima",-62.0751, Longitude)) %>%
    # Barranquillas, Colombia, geocode
    dplyr::mutate(Latitude= ifelse(Citation == "Hugot 1985" & ParasiteCorrectedName == "Trypanoxyuris clementinae" & HostCorrectedName == "Cebus albifrons", 11.00411, Latitude),
                  Longitude= ifelse(Citation == "Hugot 1985" & ParasiteCorrectedName == "Trypanoxyuris clementinae" & HostCorrectedName == "Cebus albifrons", -74.80698, Longitude)) %>%
    # Brasil, geocode
    dplyr::mutate(Latitude= ifelse(Citation == "Hugot 1985" & ParasiteCorrectedName == "Trypanoxyuris clementinae" & HostCorrectedName == "Cebus apella", -14.235, Latitude),
                  Longitude= ifelse(Citation == "Hugot 1985" & ParasiteCorrectedName == "Trypanoxyuris clementinae" & HostCorrectedName == "Cebus apella", -51.92528, Longitude)) %>%
    # ------Herder et al. 2002
    # For ungulates GMPD has Cameroon but the sampling site was Bipindi, coords updated to Bipindi, as in the primates reported in the same paper
    dplyr::mutate(Latitude= ifelse(Citation == "Herder et al. 2002"  & LocationName == "Cameroon", 3.08, Latitude),
                  Longitude= ifelse(Citation == "Herder et al. 2002"  & LocationName == "Cameroon", 10.41, Longitude)) %>%
    # ------Iori and Lanfranchi 1996
    # Somalia: the Bay, Middle and Low Scebeli regions, we geocoded the Bay region in Somalia
    dplyr::mutate(Latitude= ifelse(Citation == "Iori and Lanfranchi 1996"  & LocationName == "Somalia: the Bay, Middle and Low Scebeli regions", 2.482519, Latitude),
                  Longitude= ifelse(Citation == "Iori and Lanfranchi 1996"  & LocationName == "Somalia: the Bay, Middle and Low Scebeli regions", 43.48374, Longitude)) %>%
    #
    # Coordinates that were not assigned ecoregion, or that fall in "Rock and Ice"
    # or "Lake" categories were double-check by Paula Pappalardo looking at the
    # original citations when necessary. Most errors fixed, small islands seem
    # to have an issue but we are not including them in the analysis 
    dplyr::mutate(Latitude= ifelse(LocationName == "Yokohama, JAPAN", 35.444, Latitude),
                  Longitude= ifelse(LocationName == "Yokohama, JAPAN",  139.638, Longitude)) %>% 
    dplyr::mutate(Latitude= ifelse(LocationName == "Long Point, a peninsula originating from the north shore of Lake Erie, Ontario, Canada", 42.580, Latitude),
                  Longitude= ifelse(LocationName == "Long Point, a peninsula originating from the north shore of Lake Erie, Ontario, Canada",  -80.388, Longitude)) %>% 
    dplyr::mutate(Latitude= ifelse(LocationName == "northeastern Minnesota: parts of Cook, Lake, and Saint Louis counties in northeastern Minnesota incl", 47.494, Latitude),
                  Longitude= ifelse(LocationName == "northeastern Minnesota: parts of Cook, Lake, and Saint Louis counties in northeastern Minnesota incl",  -91.392, Longitude)) %>% # I mapped the state park Finland State forest
    dplyr::mutate(Latitude= ifelse(LocationName == "Onslow County, North Carolina", 34.65401, Latitude),
                  Longitude= ifelse(LocationName == "Onslow County, North Carolina",  -77.4702, Longitude)) %>% # I used geocode
    dplyr::mutate(Latitude= ifelse(LocationName == "Seatuck National Widlife Refuge, Town of Islip, Suffolk County, New York", 40.715, Latitude),
                  Longitude= ifelse(LocationName == "Seatuck National Widlife Refuge, Town of Islip, Suffolk County, New York",  -73.208, Longitude)) %>%
    dplyr::mutate(Latitude= ifelse(LocationName == "Los Tuxtlas", 18.484, Latitude),
                  Longitude= ifelse(LocationName == "Los Tuxtlas", -95.01, Longitude),
                  LocationName= ifelse(LocationName == "Los Tuxtlas", "Los Tuxtlas biosphere reserve", LocationName)) %>% #Trejo-Macias et al 2007, coords from Wikipedia
    dplyr::mutate(Latitude= ifelse(LocationName == "Philippines", 14.7794, Latitude),
                  Longitude= ifelse(LocationName == "Philippines",  121.4383, Longitude)) %>%  # Philippines centroid in water, fixed to Quezon area
    dplyr::mutate(Latitude= ifelse(LocationName == "Kirindy", -20.07, Latitude),
                  Longitude= ifelse(LocationName == "Kirindy", 44.6569, Longitude),
                  LocationName= ifelse(LocationName == "Kirindy", "Kirindy forest, Madagascar", LocationName)) %>% # coords from Wikipedia 
    dplyr::mutate(Latitude= ifelse(LocationName == "Yacu", 30.358611, Latitude),
                  Longitude= ifelse(LocationName == "Yacu",  130.528611, Longitude),
                  LocationName= ifelse(LocationName == "Yacu", "Yakushima Island", LocationName)) %>% # coords from Wikipedia 
    dplyr::mutate(Latitude= ifelse(LocationName == "Taiwan", 25.033333, Latitude),
                  Longitude= ifelse(LocationName == "Taiwan",  121.633333, Longitude)) %>% # coords from Wikipedia to the capitol
    dplyr::mutate(Latitude= ifelse(LocationName == "Mahambo", -17.483333, Latitude),
                  Longitude= ifelse(LocationName == "Mahambo",  49.466667, Longitude)) %>% # coords from Wikipedia
    dplyr::mutate(Latitude= ifelse(LocationName == "Lokobe", -13.399167, Latitude),
                  Longitude= ifelse(LocationName == "Lokobe",   48.318333, Longitude),
                  LocationName= ifelse(LocationName == "Lokobe", "Lokobe reserve, Madagascar", LocationName)) %>% # coords from Wikipedia
    dplyr::mutate(Latitude= ifelse(LocationName == "Antongyl, Madagascar", -15.468, Latitude),
                  Longitude= ifelse(LocationName == "Antongyl, Madagascar",  49.947, Longitude)) %>% # paper says North-east of Antongyl Bay, I assigned a point in land to the North-east of the bay
    dplyr::mutate(Latitude= ifelse(LocationName == "Itan, Taiwan", 23.69781, Latitude),
                  Longitude= ifelse(LocationName == "Itan, Taiwan",  120.960500, Longitude)) %>% # geocode coordinates for Taiwan
    dplyr::mutate(Latitude= ifelse(LocationName == "Darien Province, Panama", 8.408333, Latitude),
                  Longitude= ifelse(LocationName == "Darien Province, Panama",  -77.915, Longitude)) %>% # coords from Wikipedia
    dplyr::mutate(Latitude= ifelse(LocationName == "Kamuri", 36.34475, Latitude),
                  Longitude= ifelse(LocationName == "Kamuri",  136.937593, Longitude)) %>% # approximated in Google maps from map in paper
    dplyr::mutate(Latitude= ifelse(LocationName == "Kinkasan, Japan", 38.295278, Latitude),
                  Longitude= ifelse(LocationName == "Kinkasan, Japan",  141.566667, Longitude)) %>% # coords from Wikipedia
    dplyr::mutate(Latitude= ifelse(LocationName == "Port Dickson", 2.516667, Latitude),
                  Longitude= ifelse(LocationName == "Port Dickson",  101.8, Longitude)) %>% # coords from Wikipedia
    dplyr::mutate(Latitude= ifelse(LocationName == "Smitswinkel Bay", -34.263857, Latitude),
                  Longitude= ifelse(LocationName == "Smitswinkel Bay",  18.465880	, Longitude)) %>% # Moved a bit inland
    dplyr::mutate(Latitude= ifelse(LocationName == "Kilifi", -3.435823, Latitude),
                  Longitude= ifelse(LocationName == "Kilifi",  39.841686, Longitude)) %>% # Moved a bit inland
    dplyr::mutate(Latitude= ifelse(LocationName == "Puttalam", 8.045117, Latitude),
                  Longitude= ifelse(LocationName == "Puttalam",  79.854729, Longitude)) %>% # Moved a bit inland
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Gombe Stream National Park, Lake Tanganyika, Tanzania",  -4.666667, Latitude),
                  Longitude= ifelse(LocationName == "Gombe Stream National Park, Lake Tanganyika, Tanzania",  29.633333	, Longitude)) %>% # coords from Wikipedia
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Gibraltar",  36.138706, Latitude),
                  Longitude= ifelse(LocationName == "Gibraltar",  -5.347422	, Longitude)) %>% # coords from Wikipedia
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Sese Islands",  -0.423710, Latitude),
                  Longitude= ifelse(LocationName == "Sese Islands",  32.250590, Longitude)) %>% # mapped Bugala Island from Sese Islands
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Copenhagen",  55.567727, Latitude),
                  Longitude= ifelse(LocationName == "Copenhagen",  11.841038, Longitude)) %>%  # updated geocode
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Aomori & Akita Prefectures",  40.286904, Latitude),
                  Longitude= ifelse(LocationName == "Aomori & Akita Prefectures",  140.605264, Longitude)) %>%  # mapped a middle point between the two
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Swiss Jura formation",  47.34055, Latitude),
                  Longitude= ifelse(LocationName == "Swiss Jura formation",  7.246221, Longitude)) %>%  # updated geocode
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Assateague Island National Seashore, Worcester Co., Maryland, Accomack Co., Virginia",  38.059352, Latitude),
                  Longitude= ifelse(LocationName == "Assateague Island National Seashore, Worcester Co., Maryland, Accomack Co., Virginia",  -75.222577, Longitude)) %>%  # Moved a bit towards the middle of the island
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Yurihonjo, Gero, and Maibara, JAPAN",  35.605851, Latitude),
                  Longitude= ifelse(LocationName == "Yurihonjo, Gero, and Maibara, JAPAN",  136.753989, Longitude)) %>%  # mapped middle point between Gero and Maibara
    dplyr::mutate(Latitude= ifelse(LocationName ==  "South-central Manitoba, CANADA",  50.903889, Latitude),
                  Longitude= ifelse(LocationName == "South-central Manitoba, CANADA",  -99.513889, Longitude)) %>%  # updated coords from paper
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Grandview, Gilbert Plains, Rossburn, & Clan William Municipalities, Manitoba, CANADA",  50.831111, Latitude),
                  Longitude= ifelse(LocationName == "Grandview, Gilbert Plains, Rossburn, & Clan William Municipalities, Manitoba, CANADA",  -100.210189, Longitude)) %>%  # I mapped Riding Mountain National Park, that is in the middle of those regions
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Greenland",  66.290337, Latitude),
                  Longitude= ifelse(LocationName == "Greenland",  -51.715127, Longitude)) %>%  # I mapped a south western point to represent the area with low arctic tundra, not covered by ice, instead of the centroid for Greenland that falls in "Rock and Ice"
    dplyr::mutate(Latitude= ifelse(Citation ==  "Hersteinsson et al. 1993", 64.12652, Latitude),
                  Longitude= ifelse(Citation == "Hersteinsson et al. 1993",  -21.81744, Longitude)) %>%  # updated geocode
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Tsitsikama Forest National Park, Cape Province", -33.834701, Latitude),
                  Longitude= ifelse(LocationName == "Tsitsikama Forest National Park, Cape Province",  23.452318, Longitude)) %>%  # moved a bit inland on the Garden Route National Park, that now includes Tsitsikama Forest
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Beaufort county, South Carolina", 32.230732, Latitude),
                  Longitude= ifelse(LocationName == "Beaufort county, South Carolina",  -80.712103, Longitude)) %>%  # moved a bit inland
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Kanagawa and Fukushima prefectures, JAPAN", 35.44751, Latitude),
                  Longitude= ifelse(LocationName == "Kanagawa and Fukushima prefectures, JAPAN",  139.6423, Longitude)) %>%  # geocode for Kanagawa prefecture
    dplyr::mutate(Latitude= ifelse(LocationName ==  "St. Anthony herd, Newfoundland, Canada", 51.391652, Latitude),
                  Longitude= ifelse(LocationName == "St. Anthony herd, Newfoundland, Canada",  -55.622182, Longitude)) %>%  #  mapped based on paper map but looking in google maps for reference to not fall in water
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Bay de Verde, Newfoundland, Canada", 47.998985, Latitude),
                  Longitude= ifelse(LocationName == "Bay de Verde, Newfoundland, Canada",  -53.108995, Longitude)) %>%  #  mapped based on paper map but looking in google maps for reference to not fall in water
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Northern peninsula herd, Newfoundland, Canada", 50.730418, Latitude),
                  Longitude= ifelse(LocationName == "Northern peninsula herd, Newfoundland, Canada",  -57.035497, Longitude)) %>%  # mapped based on paper map but looking in google maps for reference to not fall in water
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Lower Columbia River, from Portland downstream to Astoria, Oregon", 45.920857, Latitude),
                  Longitude= ifelse(LocationName == "Lower Columbia River, from Portland downstream to Astoria, Oregon",  -123.373194, Longitude)) %>%  # moved inland betwen the two areas
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Kugluktuk, Nunavut, CANADA", 67.818254, Latitude),
                  Longitude= ifelse(LocationName == "Kugluktuk, Nunavut, CANADA",  -115.129287, Longitude)) %>%  # moved a bit inland 
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Cape Cross, Wolf Bay, & Atlas Bay, NAMIBIA", -27.311242, Latitude),
                  Longitude= ifelse(LocationName == "Cape Cross, Wolf Bay, & Atlas Bay, NAMIBIA",   15.562399, Longitude)) %>%  # Aproximated from the "Wolf Bay & Atlas Bay" location draw in the paper
    dplyr::mutate(Latitude= ifelse(LocationName ==  "south Andaman, Wrafters Creek, Baratang (India)",  11.925162, Latitude),
                  Longitude= ifelse(LocationName == "south Andaman, Wrafters Creek, Baratang (India)",    92.660598, Longitude)) %>%  # I mapped South Andaman island
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Kodiak Island, Alaska",  57.49125, Latitude),
                  Longitude= ifelse(LocationName == "Kodiak Island, Alaska",  -153.495, Longitude)) %>%  # updated with geocode
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Admiralty Island, Alaska",  57.68529, Latitude),
                  Longitude= ifelse(LocationName == "Admiralty Island, Alaska", -134.4895, Longitude)) %>%  # updated with geocode  
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Monhegan Island, Maine",  43.765593, Latitude),
                  Longitude= ifelse(LocationName == "Monhegan Island, Maine", -69.312324, Longitude)) %>%  # moved inland with Google Maps 
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Merritt Island, Florida", 28.343995, Latitude),
                  Longitude= ifelse(LocationName == "Merritt Island, Florida", -80.694766, Longitude)) %>%  # moved inland with Google Maps 
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Kuiu Island, Alaska", 56.56547, Latitude),
                  Longitude= ifelse(LocationName == "Kuiu Island, Alaska", -134.1334, Longitude)) %>%  # updated with geocode
    dplyr::mutate(Latitude= ifelse(LocationName ==  "San Juan Islands, Washington", 48.55137, Latitude),
                  Longitude= ifelse(LocationName == "San Juan Islands, Washington", -123.0781, Longitude)) %>%  # updated with geocode for San Juan Island instead of the archipielago
    dplyr::mutate(Latitude= ifelse(LocationName ==  "Revillagigedo Island, Alaska", 55.62712, Latitude),
                  Longitude= ifelse(LocationName == "Revillagigedo Island, Alaska", -131.4989, Longitude)) %>%  # updated with geocode 
    dplyr::mutate(Latitude= ifelse(LocationName == "Shetland, Scotland", 60.32067, Latitude),
                  Longitude= ifelse(LocationName == "Shetland, Scotland", -1.288577, Longitude)) %>%  # updated with geocode for "Mainland, Shetland Islands"
    dplyr::mutate(Latitude= ifelse(LocationName == "Hordaland, Norway", 60.27337, Latitude),
                  Longitude= ifelse(LocationName == "Hordaland, Norway", 5.722019, Longitude)) %>%  # updated with geocode 
    dplyr::mutate(Latitude= ifelse(LocationName == "Prince Edward Island", 46.324343, Latitude),
                  Longitude= ifelse(LocationName == "Prince Edward Island", -63.399925, Longitude)) %>%  # updated with geocode and moved a bit inland
    dplyr::mutate(Latitude= ifelse(LocationName == "Alaska Peninsula (Katmai coast and Black Lake)", 58.59753, Latitude),
                  Longitude= ifelse(LocationName == "Alaska Peninsula (Katmai coast and Black Lake)", -154.6937, Longitude)) %>%  # updated with geocode for Katmai National Park and Preserve
    dplyr::mutate(Latitude= ifelse(LocationName ==  "SW Greenland", 60.917365, Latitude),
                  Longitude= ifelse(LocationName ==  "SW Greenland", -46.038209, Longitude)) %>%  # They refer to SW Greenland to sites 6-8, I used google maps for site 6 "Narsaq district" geocode for 7 doesn't work
    dplyr::mutate(Latitude= ifelse(LocationName == "Port Renfrew, Vancouver Island", 48.544313, Latitude),
                  Longitude= ifelse(LocationName ==  "Port Renfrew, Vancouver Island", -124.398249, Longitude)) %>%  # updated with geocode and moved inland with google maps
    dplyr::mutate(Latitude= ifelse(LocationName == "Svalbard, Norway",  77.26259, Latitude),
                  Longitude= ifelse(LocationName ==  "Svalbard, Norway", 16.07378, Longitude)) %>%  # # updated with geocode for Sor-Spitsbergen National Park because the centroid falls in water
    dplyr::mutate(Latitude= ifelse(LocationName == "Bogesund, near Stockholm in south-central Sweden",  59.40198, Latitude),
                  Longitude= ifelse(LocationName ==  "Bogesund, near Stockholm in south-central Sweden", 18.24724, Longitude)) %>%  # updated with geocode for Bogsund, Sweden
    dplyr::mutate(Latitude= ifelse(LocationName == "Aukra, Norway",  62.819637, Latitude),
                  Longitude= ifelse(LocationName ==  "Aukra, Norway", 6.861088, Longitude)) %>%  # updated with geocode, map inland with google maps
    dplyr::mutate(Latitude= ifelse(LocationName == "Thomas Bay, Alaska",  57.046199, Latitude),
                  Longitude= ifelse(LocationName ==  "Thomas Bay, Alaska", -132.771990, Longitude)) %>%   # mapped inland
    dplyr::mutate(Latitude= ifelse(LocationName == "Munson Island, Florida",  24.62376, Latitude),
                  Longitude= ifelse(LocationName ==  "Munson Island, Florida", -81.40036, Longitude)) %>%  # updated with geocode
    dplyr::mutate(Latitude= ifelse(LocationName == "Isle au Haut, Maine",  44.044125, Latitude),
                  Longitude= ifelse(LocationName ==  "Isle au Haut, Maine", -68.621807, Longitude)) %>%  # updated with google maps
    dplyr::mutate(Latitude= ifelse(LocationName == "Vancouver Island, British Columbia",  49.65064, Latitude),
                  Longitude= ifelse(LocationName ==  "Vancouver Island, British Columbia", -125.4494, Longitude)) %>%  # updated with geocode 
    dplyr::mutate(Latitude= ifelse(LocationName == "Island of Uto, Stockholm archipelago, Sweden",  58.93612, Latitude),
                  Longitude= ifelse(LocationName ==  "Island of Uto, Stockholm archipelago, Sweden", 18.25031, Longitude)) %>%  # updated with geocode 
    dplyr::mutate(Latitude= ifelse(LocationName == "Sweden (Svealand)",  58.25279, Latitude),
                  Longitude= ifelse(LocationName ==  "Sweden (Svealand)", 13.05964, Longitude)) %>%  # updated with geocode for Svealand, Sweden                   
    dplyr::mutate(Latitude= ifelse(LocationName == "Temi River",  -6.369028, Latitude),
                  Longitude= ifelse(LocationName ==  "Temi River",34.88882, Longitude)) %>%  # northern Tanzania according to paper, geocode to Tanzania 
    dplyr::mutate(Latitude= ifelse(LocationName == "Geruma Island, Japan",   26.18171, Latitude),
                  Longitude= ifelse(LocationName ==  "Geruma Island, Japan", 127.2893, Longitude)) %>%  # updated with geocode
    dplyr::mutate(Latitude= ifelse(LocationName == "Rogaland county, Norway",  59.210387, Latitude),
                  Longitude= ifelse(LocationName ==  "Rogaland county, Norway", 6.221349, Longitude)) %>%  # geocode falls in water, mapped inland with google maps
    dplyr::mutate(Latitude= ifelse(LocationName == "Svanoy, Norway",   61.49005, Latitude),
                  Longitude= ifelse(LocationName ==  "Svanoy, Norway", 5.145064, Longitude)) %>% # geocode falls in water,  mapped inland with google maps
    dplyr::mutate(Latitude= ifelse(LocationName == "Vega, Norway",   65.644373, Latitude),
                  Longitude= ifelse(LocationName ==  "Vega, Norway",  11.887690, Longitude)) %>% # geocode falls in water, mapped in center of Vega Island with google maps
    dplyr::mutate(Latitude= ifelse(LocationName == "Portal, Spain",   43.462060, Latitude),
                  Longitude= ifelse(LocationName ==  "Portal, Spain",  -5.647646, Longitude)) %>% # approximated in google maps from map in paper
    dplyr::mutate(Latitude= ifelse(LocationName == "Eaglemont near Buffalo Range Ranch",  -21.0167, Latitude),
                  Longitude= ifelse(LocationName ==  "Eaglemont near Buffalo Range Ranch", 31.5833, Longitude)) %>% # Coordinates from the web (https://travelingluck.com/Africa/Zimbabwe/Zimbabwe+(general)/_6297106_Buffalo+Range.html) for Buffalo Range airport, checked with google maps and map drawed in paper
    dplyr::mutate(Latitude= ifelse(LocationName == "Eastern NORWAY",   59.56844, Latitude),
                  Longitude= ifelse(LocationName ==  "Eastern NORWAY",  9.279733, Longitude))  %>%# Eastern NORWAY point falls in Australia, corrected to geocode (ggmap) for "eastern Norway"
    dplyr::mutate(Latitude= ifelse(LocationName == "Lake Albert"|LocationName == "Kawa, Lake Albert, Africa",   2.119025, Latitude),
                  Longitude= ifelse(LocationName ==  "Lake Albert"|LocationName == "Kawa, Lake Albert, Africa",  31.526991, Longitude)) %>% # moved a bit inland
    dplyr::mutate(Latitude= ifelse(LocationName == "Nemuro, Hokkaido"|LocationName == "Nemuro district, Hokkaido Japan",   43.231560, Latitude),
                  Longitude= ifelse(LocationName ==  "Nemuro, Hokkaido"|LocationName == "Nemuro district, Hokkaido Japan",   145.478401, Longitude)) %>% # I mapped a point inland in Nemuro District
    dplyr::mutate(Latitude= ifelse(LocationName == "Sweden"|LocationName == "south and central Sweden",   60.112100, Latitude),
                  Longitude= ifelse(LocationName ==  "Sweden"|LocationName == "south and central Sweden",   18.622301, Longitude))  %>%# updated geocode, moved inland because it still fails
    dplyr::mutate(Latitude= ifelse(LocationName == "Avalon herd, Newfoundland, Canada"|LocationName == "Avalon, Newfoundland, Canada",   46.910939, Latitude),
                  Longitude= ifelse(LocationName ==  "Avalon herd, Newfoundland, Canada"|LocationName == "Avalon, Newfoundland, Canada",   -53.308928, Longitude))  %>% # mapped based on paper map, but looking in google maps for reference to not fall in water
    dplyr::mutate(Latitude= ifelse(LocationName == "Altiplano region, north of Murcia, Spain"|Citation == "Martinez-Carrasco et al. 2007",    38.25, Latitude),
                  Longitude= ifelse(LocationName ==  "Altiplano region, north of Murcia, Spain"|Citation == "Martinez-Carrasco et al. 2007",   -1.045, Longitude))  %>% # found by Claire
    dplyr::mutate(Latitude= ifelse(LocationName == "West Coast National Park"|Citation =="Boomker et al 2000",    -33.17, Latitude),
                  Longitude= ifelse(LocationName == "West Coast National Park"|Citation == "Boomker et al 2000",   18.149, Longitude))  %>% # found by Claire
    # Bentler et al. 2007 double-checked, it mentioned two specific locations on Erie county the 
    # U.S. National Aeronautics and Space Administration's Plum Brook Station (PBS)
    # Old Woman Creek National Estuarine Research Reserve (OWC)(we used this one)
    dplyr::mutate(Latitude= ifelse(LocationName == "Erie county, Ohio"|Citation == "Bentler et al. 2007",    41.377 , Latitude),
                  Longitude= ifelse(LocationName ==  "Erie county, Ohio"|Citation == "Bentler et al. 2007",   -82.509, Longitude))  %>% # 
    dplyr::mutate(Latitude= ifelse(LocationName == "Maje, Panama"|Citation == "Seymour et al. 1983",     9.2536, Latitude),
                  Longitude= ifelse(LocationName == "Maje, Panama"|Citation ==  "Seymour et al. 1983",  -78.9341, Longitude))  %>%  # point approximated from map in paper
    # longitude had the wrong sign, it should be west; coords from paper still fall in water
    # changed to wikipedia coords for Santa Marta == Los Tuxclas, Mexico
    dplyr::mutate(Latitude= ifelse(Citation == "Valdespino et al. 2010"|Citation == "Rodhain 1941",     9.2536, Latitude),
                  Longitude= ifelse(Citation == "Valdespino et al. 2010"|Citation ==  "Rodhain 1941",  -78.9341, Longitude))  %>%
    # external info from Rev. Inst. Med. trop. S. Paulo vol. 39 no. 2
    #São Paulo Mar./Apr. 1997 http://dx.doi.org/10.1590/S0036-46651997000200011 
    #By Guillermo M. DENEGRI (1) & Jorge PEREZ-SERRANO (2)
    # New location is  "San Cayetano, Corrientes, Argentina"
    dplyr::mutate(Latitude= ifelse(LocationName == "Argentina"|Citation == "Stuart et al. 1998",     -27.49462, Latitude),
                  Longitude= ifelse(LocationName == "Argentina"|Citation ==  "Stuart et al. 1998",  -58.79827, Longitude),
                  LocationName= ifelse(LocationName == "Argentina"|Citation == "Stuart et al. 1998", "San Cayetano, Corrientes, Argentina", LocationName))   # 
  
  # return dataset with all the coordinates fixed
  return(edited.df)
}
#
limitOccsToIUCNrange <- function(occs, polys, buffer){
  # Keeps only parasite occurrences within a host IUCN range
  #
  # Args:
  #   occs: dataframe with columns "Longitude" and "Latitude"
  #   gmpd.mammals: shape file with IUCN ranges for hosts in GMPD
  #   buffer: distance in km to buffer the points
  # Returns:
  #   Same dataframe with IUCN columns added, where rows outside of range have been filtered out
  # 
  library(rgeos)
  parpoints <- SpatialPoints(coords= data.frame(lon= occs$Longitude, lat= occs$Latitude),
                             proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  # print host name for debugging
  #print(occs$HostCorrectedName[1])
  # subset hosts that are in our database
  polys.in.DB <- which(as.character(polys@data$binomial) %in% occs$HostCorrectedName[1])
  this.host <- polys[polys.in.DB,]
  # transform points and host ranges to equal area projection
  eaproj <- CRS("+proj=laea +lat_0=0 +lon_0=0 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  parpoints.ea <- spTransform(parpoints, eaproj)
  host.ea <- spTransform(this.host, eaproj)
  # create a buffer around the points
  bufferwidth <- buffer * 1000
  parpoints.buffer <- gBuffer(parpoints.ea, byid = T , width = bufferwidth) 
  # overlay points on host IUCN range
  iucn.points <- over(x= parpoints.buffer, y= host.ea)
  # merge with points info 
  iucn.occs <- cbind(occs, iucn.points)
  # check how it works
  #tosee <- spTransform(parpoints.buffer, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  # remove rows without IUCN information
  iucn.occs <- iucn.occs %>% filter(!is.na(binomial))
  # return dataframe
  return(iucn.occs)
}
#
matchGMPDwithIUCN <- function(mydf){
  # Replace host names to match IUCN names
  #
  # Args:
  #   mydf: dataframe with columns "HostCorrectedName"
  # Returns:
  #   same dataframe with issues fixed
  new.df <- mydf %>%   # 
    # carnivores and ungulates 
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Alcelaphus lichtensteinii", "Alcelaphus buselaphus", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Bos frontalis", "Bos gaurus", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Bos grunniens", "Bos mutus", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Capra hircus", "Capra aegagrus", HostCorrectedName)) %>%
    #dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Equus asinus", "Equus africanus", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Equus burchellii", "Equus quagga", HostCorrectedName)) %>%
    #dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Equus caballus", "Equus ferus", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Felis manul", "Otocolobus manul", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Lama glama", "Lama guanicoe", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Leopardus pajeros", "Leopardus colocolo", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Neotragus moschatus", "Nesotragus moschatus", HostCorrectedName)) %>%
    #dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Ovis aries", "Ovis orientalis", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Puma yagouaroundi", "Herpailurus yagouaroundi", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Taurotragus oryx", "Tragelaphus oryx", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Alces americanus", "Alces alces", HostCorrectedName)) %>%
    # primates 
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Callithrix emiliae", "Mico emiliae", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Callithrix argentata", "Mico argentatus", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Callithrix pygmaea", "Cebuella pygmaea", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Cebus apella",  "Sapajus apella", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Cebus libidinosus",  "Sapajus libidinosus", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Cebus xanthosternos",  "Sapajus xanthosternos", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Cercopithecus albogularis",  "Cercopithecus mitis", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Cercopithecus kandti",  "Cercopithecus mitis", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Galago demidoff",  "Galagoides demidoff", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Galago zanzibaricus",  "Galagoides zanzibaricus", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Cercopithecus lhoesti",  "Allochrocebus lhoesti", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Cercopithecus lhoesti",  "Allochrocebus lhoesti", HostCorrectedName)) %>%
    dplyr::mutate(HostCorrectedName= ifelse(HostCorrectedName == "Mico chrysoleucus",  "Mico chrysoleucos", HostCorrectedName))
  return(new.df)  
}


# ------Random useful functions
# -----------------------------
# Functions used in different steps of data analyis
#
coords2country <- function(mydf){  
  # Assigns country to geographic occurrences
  # (most of the code here was found in online examples)
  #
  # Args:
  #   dataframe with columns "Longitude" and "Latitude"
  # Returns:
  #   The same dataframe with an extra colum named "Country"
  library(sp)
  library(rworldmap)
  myres <- mydf
  # get world map
  countriesSP <- getMap(resolution='low')
  #countriesSP <- getMap(resolution='high') #you could use high res map from rworldxtra if you were concerned about detail
  # create a spatial points data frame with same projection that the maps we use to query the points
  coordinates(mydf) <- ~ Longitude + Latitude
  proj4string(mydf) <- proj4string(countriesSP)
  # convert our list of points to a SpatialPoints object
  # use 'over' to get indices of the Polygons object containing each point 
  indices <- over(mydf, countriesSP)
  # return the ADMIN names of each country
  myres$country <- indices$ADMIN  
  #indices$ISO3 # returns the ISO3 code 
  myres$continent <- indices$continent   # returns the continent (6 continent model)
  myres$region <- indices$REGION   # returns the continent (7 continent model)
  return(myres)
}
#
selectUSoccs <- function(occs){
  # Creates a dataset of parasites occurrences only in the continental United States
  #
  # Args:
  #   mydf: cleaned occurrences file
  # Returns:
  #   Dataframe with occurrences only in continental United States
  library(dplyr)
  # subset ungulates and carnivores
  unca <- occs %>% dplyr::filter(HostOrder == "Artiodactyla"|HostOrder == "Carnivora"|HostOrder == "Perissodactyla") %>% 
                   coords2country() %>% 
                   dplyr::filter(country == "United States of America") %>% 
                   dplyr::filter(Latitude < 53) %>%
                   dplyr::filter(Longitude > -140)
}
#
selectAFRICAoccs <- function(occs){
  # Creates a dataset of parasites occurrences only in the continental United States
  #
  # Args:
  #   mydf: cleaned occurrences file
  # Returns:
  #   Dataframe with occurrences only in continental United States
  library(dplyr)
  # subset ungulates and carnivores
  pri <- occs %>% dplyr::filter(HostOrder == "Primates") %>% 
    coords2country() %>% 
    dplyr::filter(continent == "Africa") 
}
#
runNN2SearchOnGrid <- function(mydf, exp.grid, n.neighbours, radius){
  # Assigns ID to occurrence points on a grid using the function "nn2" (neighbour joining algorithm).
  #
  # Args:
  #   mydf= dataframe with parasite coordinates (Longitude and Latitude columns) and parasite names (ParasiteCorrectedName column)
  #   exp.grid= expanded grid to run the search on (lon, lat columns)
  #   n.neighbours= number of neighbour cells to include 
  #   radius= distance to search for neighbors, in decimal degrees
  # Returns:
  #   A dataframe with the parasite occurrences and a column for their ids on the grid 
  #
  # Load libraries:
  library(RANN)
  library(dplyr)
  library(magrittr)
  #
  # --------NEAREST NEIGHBOUR SEARCH-------------
  # find the closest point on grid for each of our occurrence points
  # run a nearest neighbour search with function "nn2", the result is a list with the first component
  # being the id in the grid dataframe and the second component being the distance to the grid
  # assign id to coordinates to merge later
  mydf$coord.id <- seq(1:nrow(mydf))
  # coordinates to interpolate
  parcoords <- mydf[,c("Longitude", "Latitude")]
  # create list to hold id results
  idlist <- list()
  # loop through each coordinate to get neighbour cells ids
  for (i in 1:nrow(parcoords)){
    idname <- paste(i)
    print(paste("working in row ", i))
    # run neighbor joining search 
    closepoints <- nn2(query= parcoords[i,], data= exp.grid[,c("lon", "lat")], k= n.neighbours, radius= radius, searchtype=c("radius"))
    # identify point to merge the closest cells id 
    id.df <- data.frame(coord.id= rep(i, n.neighbours), ids= as.numeric(closepoints[[1]]), dist= as.numeric(closepoints[[2]]))
    id.df <- id.df %>% filter(ids!=0)
    idlist[[idname]] <- id.df
  }
  # combine all id dataframes
  fullids <- do.call("rbind", idlist)
  # merge with the parasite occurrences
  mydf.1 <- join(mydf, fullids, by= "coord.id", type= "right")
  return(mydf.1)
}
#
occs2ecoregions_EqualArea <- function(occs, ecos){  
  # Assigns ecoregion to geographic occurrences
  #
  # Args:
  #   occs: dataframe with columns "Longitude" and "Latitude"
  #   ecos: shape file with terrestrial ecoregions
  # Returns:
  #   A dataframe with ecoregion name (ECO_NAME), parcounts, hostcounts,
  #   citation counts and occurrences counts in each ecoregion
  #
  # load libraries
  library(raster)
  library(rgdal)
  library(plyr)
  library(ggplot2)
  library(GGally)
  
  # -----Overlay occurrences and ecoregions-----------------
  # create Spatial Points objects with the parasite coordinates, assigning exactly the same projection that the ecoregions have.
  parpoints <- SpatialPoints(coords= data.frame(lon= occs$Longitude, lat= occs$Latitude), proj4string = CRS(ecos@proj4string@projargs))
  parpoints.ea <- spTransform(parpoints, CRS("+proj=laea +lat_0=0 +lon_0=0 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  # transform ecoregions object
  ecos.ea <- spTransform(ecos, CRS("+proj=laea +lat_0=0 +lon_0=0 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  # overlay parasite occurrences over ecoregions: x = "SpatialPoints", y = "SpatialPolygonsDataFrame". Returns a data.frame of the second argument with row entries corresponding to the first argument
  ecopoints <- over(x = parpoints.ea, y = ecos.ea)
  #nrow(ecopoints[is.na(ecopoints$ECO_ID),]) # check how many NA (a ecoregion was not assigned) 
  # keep ecoregions useful columns
  ecopoints.1 <- subset(ecopoints, select= c(ECO_NAME, ECO_ID, SHAPE_AREA))
  # add ecoregion information to parasite occurrences dataframe
  occs.1 <- cbind(occs, ecopoints.1)
  # save temporary object occs.1 to use in accumulation curves analysis
  #filename <- paste(here::here(),"/Files_output/datasets/OccsByEcoregion_world_", mg, ".csv", sep="")
  write.csv(occs.1, filename, row.names= F)
  # ----Extra check--------
  # There are two categories in the ecoregions dataframe that are not a proper ecoregions: "Rock and Ice" and "Lake". 
  #View(occs.1[which(occs.1$ECO_NAME == "Rock and Ice"),])
  #View(occs.1[which(occs.1$ECO_NAME == "Lake"),]) 
  # Filter these two categories 
  occs.2 <- occs.1 %>% dplyr::filter(!(ECO_NAME == "Rock and Ice"|ECO_NAME == "Lake")) %>%
                       dplyr::filter(!is.na(ECO_ID))

  # ---------Parasite Species Richness-------------
  # count parasite richness and other info in each grid cell
  summbyeco <- ddply(occs.2, .(factor(ECO_ID), ECO_NAME, SHAPE_AREA),summarize, parcounts = length(unique(as.character(ParasiteCorrectedName))), hostcounts = length(unique(as.character(HostCorrectedName))), citcounts = length(unique(as.character(Citation))), occscounts = length(Latitude))
  # subset info to merge later
  pardata <- subset(summbyeco, select= c(ECO_NAME, parcounts, hostcounts, citcounts, occscounts))
  pardata$ECO_NAME <- factor(pardata$ECO_NAME)
  return(pardata)
  
  # remove things not in use
  rm(ecopoints, ecopoints.1, occs, occs.1, occs.2, nas, noteco, parpoints, pardata, summbyeco)
}
#
occs2ecoregions_EqualArea_US <- function(occs, ecos){  
  # Assigns ecoregion to geographic occurrences
  # US centroid lat= 39.828175, lon= -98.5795
  # Args:
  #   occs: dataframe with columns "Longitude" and "Latitude"
  #   ecos: shape file with terrestrial ecoregions
  # Returns:
  #   A dataframe with ecoregion name (ECO_NAME), parcounts, hostcounts,
  #   citation counts and occurrences counts in each ecoregion
  #
  # load libraries
  library(raster)
  library(rgdal)
  library(plyr)
  library(ggplot2)
  library(GGally)
  
  # -----Overlay occurrences and ecoregions-----------------
  # create Spatial Points objects with the parasite coordinates, assigning exactly the same projection that the ecoregions have.
  parpoints <- SpatialPoints(coords= data.frame(lon= occs$Longitude, lat= occs$Latitude), proj4string = CRS(ecos@proj4string@projargs))
  # transform points to equal area projection
  parpoints.ea <- spTransform(parpoints, CRS("+proj=laea +lat_0=39.828175 +lon_0=-98.5795 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  # transform ecoregions object to equal area projection (it takes a little while)
  ecos.ea <- spTransform(ecos, CRS("+proj=laea +lat_0=39.828175 +lon_0=-98.5795 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  # overlay parasite occurrences over ecoregions: x = "SpatialPoints", y = "SpatialPolygonsDataFrame". Returns a data.frame of the second argument with row entries corresponding to the first argument
  ecopoints <- over(x = parpoints.ea, y = ecos.ea)
  nrow(ecopoints[is.na(ecopoints$ECO_ID),]) # check how many NA (a ecoregion was not assigned)
  # keep ecoregions useful columns
  ecopoints.1 <- subset(ecopoints, select= c(ECO_NAME, ECO_ID, SHAPE_AREA))
  # add ecoregion information to parasite occurrences dataframe
  occs.1 <- cbind(occs, ecopoints.1)
  # save temporary object occs.1 to use in accumulation curves analysis
  #filename <- paste(here::here(),"/Files_output/datasets/OccsByEcoregion_unca_10occs_US.csv", sep="")
  #write.csv(occs.1, filename, row.names= F)
  # ----No "bad" points in the US dataset--------
  #nrow(occs.1[which(occs.1$ECO_NAME == "Rock and Ice"),])
  #nrow(occs.1[which(occs.1$ECO_NAME == "Lake"),])

  # ---------Parasite Species Richness-------------
  # count parasite richness and other info in each grid cell
  summbyeco <- ddply(occs.1, .(factor(ECO_ID), ECO_NAME, SHAPE_AREA),summarize, parcounts = length(unique(as.character(ParasiteCorrectedName))), hostcounts = length(unique(as.character(HostCorrectedName))), citcounts = length(unique(as.character(Citation))), occscounts = length(Latitude))
  # subset info to merge later
  pardata <- subset(summbyeco, select= c(ECO_NAME, SHAPE_AREA, parcounts, hostcounts, citcounts, occscounts))
  pardata$ECO_NAME <- factor(pardata$ECO_NAME)
  return(pardata)
  
  # remove things not in use
  rm(ecopoints, ecopoints.1, occs, occs.1,  nas, noteco, parpoints, pardata, summbyeco)
}
#
occs2ecoregions_EqualArea_AF <- function(occs, ecos){  
  # Assigns ecoregion to geographic occurrences
  # central Africa centroid aprox lat= 0, lon= 24
  # Args:
  #   occs: dataframe with columns "Longitude" and "Latitude"
  #   ecos: shape file with terrestrial ecoregions
  # Returns:
  #   A dataframe with ecoregion name (ECO_NAME), parcounts, hostcounts,
  #   citation counts and occurrences counts in each ecoregion
  #
  # load libraries
  library(raster)
  library(rgdal)
  library(plyr)
  library(ggplot2)
  library(GGally)
  
  # -----Overlay occurrences and ecoregions-----------------
  # create Spatial Points objects with the parasite coordinates, assigning exactly the same projection that the ecoregions have.
  parpoints <- SpatialPoints(coords= data.frame(lon= occs$Longitude, lat= occs$Latitude), proj4string = CRS(ecos@proj4string@projargs))
  # transform points to equal area projection
  parpoints.ea <- spTransform(parpoints, CRS("+proj=laea +lat_0=0 +lon_0=24 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  # transform ecoregions object to equal area projection (it takes a little while)
  ecos.ea <- spTransform(ecos, CRS("+proj=laea +lat_0=0 +lon_0=24 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  # overlay parasite occurrences over ecoregions: x = "SpatialPoints", y = "SpatialPolygonsDataFrame". Returns a data.frame of the second argument with row entries corresponding to the first argument
  ecopoints <- over(x = parpoints.ea, y = ecos.ea)
  nrow(ecopoints[is.na(ecopoints$ECO_ID),]) # check how many NA (a ecoregion was not assigned)
  # keep ecoregions useful columns
  ecopoints.1 <- subset(ecopoints, select= c(ECO_NAME, ECO_ID, SHAPE_AREA))
  # add ecoregion information to parasite occurrences dataframe
  occs.1 <- cbind(occs, ecopoints.1)
  # save temporary object occs.1 to use in accumulation curves analysis
  #filename <- paste(here::here(),"/Files_output/datasets/OccsByEcoregion_primates_iucn50_Africa.csv", sep="")
  #write.csv(occs.1, filename, row.names= F)

  # ---------Parasite Species Richness-------------
  # count parasite richness and other info in each grid cell
  summbyeco <- ddply(occs.1, .(factor(ECO_ID), ECO_NAME, SHAPE_AREA),summarize, parcounts = length(unique(as.character(ParasiteCorrectedName))), hostcounts = length(unique(as.character(HostCorrectedName))), citcounts = length(unique(as.character(Citation))), occscounts = length(Latitude))
  # subset info to merge later
  pardata <- subset(summbyeco, select= c(ECO_NAME, SHAPE_AREA, parcounts, hostcounts, citcounts, occscounts))
  pardata$ECO_NAME <- factor(pardata$ECO_NAME)
  return(pardata)
  
  # remove things not in use
  rm(ecopoints, ecopoints.1, occs, occs.1,  nas, noteco, parpoints, pardata, summbyeco)
}
#
getParDensityLayer_US <- function (mydata, parname, bandwidth, threshold, gs, usshape){
  # Function to perform density estimation and get raster layer from it
  #
  # Args:
  # mydata: dataframe with one or more parasite occurrences in columns Longitude and Latitude.
  # parname: a character indicating the parasite name to subset
  # bandwith: bandwith of the kernel for the density estimation, is an integer
  # treshold: a quantile to use as treshold (such as 0.50, 0.75)
  # gs: grid size that we want
  # usshape: SpatialPolygonsDataFrame file with continental US
  #  
  # Returns:
  # A raster layer with the presence/absence of parasites based on density.  
  # load libraries
  library(MASS)
  library(raster)
  library(sp)
  library(rgdal)
  library(KernSmooth)
  
  # subset dataframe for only this parasite  
  mydf <- subset(mydata, mydata$ParasiteCorrectedName == parname)
  # create SpatialPoints object
  mypoints <- SpatialPoints(coords= data.frame(lon= mydf$Longitude, lat= mydf$Latitude), proj4string= CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +no_defs"))
  # create vector to represent bandwith
  mybw <- c(bandwidth, bandwidth)
  # compute the 2D binned kernel density estimate
  denkernel <- bkde2D(coordinates(mypoints), 
                      bandwidth= mybw, 
                      gridsize= c(65/gs, 35/gs),
                      range.x= list(c(-130,-65), c(20, 55)))
  denkernel$fhat[denkernel$fhat < 0.00001] <- 0 # ignore very small values
  # create raster
  denraster = raster(list(x=denkernel$x1, y=denkernel$x2, z=denkernel$fhat))
  projection(denraster) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +no_defs ")
  xmin(denraster) <- -130
  xmax(denraster) <- -65
  ymin(denraster) <- 20
  ymax(denraster) <- 55
  # assign zero density to NA
  denraster[denraster == 0] <- NA
  # define threshold
  th <- quantile(values(denraster)[!is.na(values(denraster))], threshold)
  # replace values lower than threshold with NA
  denraster[denraster < th] <- NA
  # replace values higher than threshold with 1
  denraster[denraster >= th] <- 1
  # mask with us shape to avoid density points in water
  denrasterm <- mask(denraster, usshape)
  # return raster with presences for this parasite
  return(denrasterm)
}
#
getParDensityLayer_Africa <- function (mydata, parname, bandwidth, threshold, gs, africashape){
  # Function to perform density estimation and get raster layer from it
  #
  # Args:
  # mydata: dataframe with one or more parasite occurrences in columns Longitude and Latitude.
  # parname: a character indicating the parasite name to subset
  # bandwith: bandwith of the kernel for the density estimation, is an integer
  # treshold: a quantile to use as treshold (such as 0.50, 0.75)
  # gs: grid size that we want
  # usshape: SpatialPolygonsDataFrame file with continental US
  #  
  # Returns:
  # A raster layer with the presence/absence of parasites based on density.  
  # load libraries
  library(MASS)
  library(raster)
  library(sp)
  library(rgdal)
  library(KernSmooth)
  
  # subset dataframe for only this parasite  
  mydf <- subset(mydata, mydata$ParasiteCorrectedName == parname)
  # create SpatialPoints object
  mypoints <- SpatialPoints(coords= data.frame(lon= mydf$Longitude, lat= mydf$Latitude), proj4string= CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +no_defs"))
  # create vector to represent bandwith
  mybw <- c(bandwidth, bandwidth)
  # compute the 2D binned kernel density estimate
  denkernel <- bkde2D(coordinates(mypoints), 
                      bandwidth= mybw, 
                      gridsize= c(90/gs, 80/gs),
                      range.x= list(c(-30, 60), c(-40, 40)))

  denkernel$fhat[denkernel$fhat < 0.00001] <- 0 # ignore very small values
  # create raster
  denraster = raster(list(x=denkernel$x1, y=denkernel$x2, z=denkernel$fhat))
  projection(denraster) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +no_defs ")
  xmin(denraster) <- -30
  xmax(denraster) <- 60
  ymin(denraster) <- -40
  ymax(denraster) <- 40
  # assign zero density to NA
  denraster[denraster == 0] <- NA
  # define threshold
  th <- quantile(values(denraster)[!is.na(values(denraster))], threshold)
  # replace values lower than threshold with NA
  denraster[denraster < th] <- NA
  # replace values higher than threshold with 1
  denraster[denraster >= th] <- 1
  # mask with us shape to avoid density points in water
  denrasterm <- mask(denraster, africashape)
  # return raster with presences for this parasite
  return(denrasterm)
}
#
getParDensityLayer <- function (group, parname, bandwidth, threshold, gs){
  # Function to perform density estimation and get raster layer from it
  #
  # Args:
  # group: dataframe with one or more parasite occurrences in columns Longitude and Latitude.
  # parname: a character indicating the parasite name to subset
  # bandwith: bandwith of the kernel for the density estimation, is an integer
  # treshold: a quantile to use as treshold (such as 0.50, 0.75)
  # gs: grid size that we want
  #  
  # Returns:
  # A raster layer with the presence/absence of parasites based on density.  
  # load libraries
  library(MASS)
  library(raster)
  library(sp)
  library(maptools)
  library(rgdal)
  library("KernSmooth")
  # load world map
  data("wrld_simpl")
  # Worldmap (without Antarctica)
  myworld <- wrld_simpl[wrld_simpl@data$UN!="10",]
  # subset dataframe for only this parasite  
  mydf <- subset(group, group$ParasiteCorrectedName == parname)
  # create SpatialPoints object
  mypoints <- SpatialPoints(coords= data.frame(lon= mydf$Longitude, lat= mydf$Latitude), proj4string= CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +no_defs"))
  # create vector to represent bandwith
  mybw <- c(bandwidth, bandwidth)
  # compute the 2D binned kernel density estimate
  denkernel <- bkde2D(coordinates(mypoints), 
                      bandwidth= mybw, # we need to justify this choice...
                      gridsize= c(360/gs, 180/gs),
                      range.x= list(c(-180,180), c(-90,90)))
  # test how it looks
  #contour(denkernel$x1, denkernel$x2, denkernel$fhat)
  denkernel$fhat[denkernel$fhat < 0.00001] <- 0 # ignore very small values
  # create raster
  denraster = raster(list(x=denkernel$x1, y=denkernel$x2, z=denkernel$fhat))
  projection(denraster) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +no_defs ")
  xmin(denraster) <- -180
  xmax(denraster) <- 180
  ymin(denraster) <- -90
  ymax(denraster) <- 90
  # assign zero density to NA
  denraster[denraster == 0] <- NA
  # define threshold
  th <- quantile(values(denraster)[!is.na(values(denraster))], threshold)
  # replace values lower than threshold with NA
  denraster[denraster < th] <- NA
  # replace values higher than threshold with 1
  denraster[denraster >= th] <- 1
  # mask with world to avoid density points in water
  denrasterm <- mask(denraster, myworld)
  # return raster with presences for this parasite
  return(denrasterm)
}

# ------Functions to calculate parasite species richness
# ------------------------------------------------------
#
getRichnessFromHost_Stack <- function(mydata, mamm.sub, r){
  # Calculates parasite species richness from IUCN host ranges
  #
  # Args:
  #   mydf: dataframe with columns "HostCorrectedName" and "ParasiteCorrectedName"
  #   mamm.sub: subset of IUCN polygons that correspond to gmpd host
  #   r: basic raster to use
  #   
  # Returns:
  #   raster file with parasite species richness
  #
  library(sf)
  library(sp)
  library(fasterize)
  library(dplyr)
  
  # --------------get Par richness by host filling
  mypars <- sort(unique(mydata$ParasiteCorrectedName)) # make a list of parasites
  # go one by one for each parasite aggregating the range of its hosts
  x <- stack()
  for(i in mypars){
    thispar <- mydata %>% dplyr::filter(ParasiteCorrectedName == i)
    # tag parasite for debugging 
    print(paste("working on parasite:", i))
    hosts <- unique(thispar$HostCorrectedName)
    # subset hosts
    id <- which(as.character(mamm.sub@data$binomial) %in% hosts)
    polys <- mamm.sub[id,]
    st.polys <- st_as_sf(polys)
    r.mamm <- fasterize(st.polys, r, fun="sum")
    r.mamm[values(r.mamm) > 1] <- 1 # this gives a warning when the host range
    # is not in the US, ignore, it does nothing
    # add layer for host i to a raster stack
    x <- stack(x, r.mamm)
  }
  
  # combine all raster layers into one by adding them
  r.hosts <- stackApply(x, indices= rep(1, length(mypars)), fun= sum)
  # set zeros to NAs
  r.hosts[r.hosts == 0] <- NA 
  # return raster with parasite richness
  return(r.hosts)
}
#
createDensityRasters_world <- function(mydata, mybw, myth, mg, tomask, gs){
  # Create species richness rasters with density method 
  # Args:
  #   mydata: dataframe with taxonomic group to analyze
  #   mybw: bandwidth for kernel density
  #   myth: quantile of density to set likelihood of presence and define layer to presence/absence 
  #   gridsize: grid size
  #   mg: name of taxonomic work for files
  #   group: dataframe of the taxonomic group selected for analysis
  #   tomask: SpatialPolygonDataFrame to mask layers
  #   
  # Returns:
  #   Raster layer from density method 
  #
  # get parasite names
  mypars <- sort(unique(as.character(mydata$ParasiteCorrectedName)))
  # first create an empty stack object
  x <- stack()
  # run loop to calculate density layer for each parasite and put each layer in a raster stack
  for(parname in mypars){
    # working on parasite 
    print(paste("working on parasite:", parname))
    # get a raster layer from the 2D density estimation
    myraster <- getParDensityLayer(group= mydata, parname= parname, bandwidth = mybw, threshold= myth, gs=gs)
    # add layer for parasite i to a raster stack
    x <- stack(x, myraster)
  }
  rm(myraster)
  # plot(x) # check that layers are different and look ok
  # combine all raster layers into one by adding them
  r.density.stack <- stackApply(x, indices= rep(1, length(mypars)), fun= sum)
  # assign zero density to NA (the stackApply returns zeros)
  r.density.stack[r.density.stack == 0] <- NA
  #mask layer with the world to avoid cells in water
  world.density <- mask(r.density.stack, tomask)
  # create filename for layer
  rastername <- paste(here::here(),"/Files_output/layers/world/ByDensity_", gs, "deg_", mg,"_",mybw,"bw_", myth,"th.Rdata", sep="")
  ## save final raster from grid to folder
  save(world.density, file= rastername)
  return(world.density)
}


#-----------------------------------------------------------------------------------
#--------------------------ACCUM CURVES----------------------------------------------
# function to extract data to build rarefaction curves for host parasites from GMPD
# Written by Ignacio Morales-Castilla

#' For specified parameters, subsamples host-parasites occurrences and extracts data to build accum. curves
#'
#' @param dat data.frame containing host-parasite data - e.g. GMPD
#' @param nrun scalar specifying number of iterations to subsample for each number of units (hosts or coordinates)
#' @param by.host default = 'hosts', either 'hosts' or 'coords', specifies the sampling unit
#' @param extent default = 'all', either 'all' or 'limited', sampling everything?
#' @param lonmax;lonmin;latmax;latmin specifies the extent if extent == "limited"
#'
#' @returns data.frame of dimensions n.units x nrun, to draw rarefaction curve with uncertainty

rare.curve.prop <- function(dat,nrun,by.host = c('hosts','coords','ecoreg'),extent=c("all","limited"),
                          lonmax = NULL,lonmin = NULL,latmax = NULL,latmin = NULL){
  
  ## extract unique records of hosts, parasites and coordinates
  list.hosts <- unique(dat$HostCorrectedName)
  list.parasites <- unique(dat$ParasiteCorrectedName)
  list.coords <- unique(paste(dat$Longitude,dat$Latitude))
  nhosts <- length(list.hosts)
  nparasites <- length(list.parasites)
  ncoords <- length(list.coords)
  ncoords <- round(ncoords/2,0)
  list.ecoregs <- dat$ECO_ID
  necoregs <- length(unique(list.ecoregs))
  
  
  ## subset dataset according to extent, when extent is limited only
  if(extent == 'limited'){
    dat <- subset(dat,Longitude>lonmin & Longitude<lonmax 
                & Latitude>latmin & Latitude<latmax)
    list.hosts <- unique(dat$HostCorrectedName)
    list.parasites <- unique(dat$ParasiteCorrectedName)
    list.coords <- unique(paste(dat$Longitude,dat$Latitude))
    list.ecoregs <- sort(unique(dat$ECO_ID))
    nhosts <- length(list.hosts)
    nparasites <- length(list.parasites)
    ncoords <- length(list.coords)
    ncoords <- round(ncoords/2,0)
    list.ecoregs <- dat$ECO_ID
    necoregs <- length(unique(list.ecoregs))
  }
  
  # method to count parasites within a subset of host species
  if(by.host == 'hosts'){
    
    # generate array to store accumulation curves 
    rare.data <- array(NA,dim = c(nhosts,nrun+1))
    
    # loop along hosts
    for(i in 1:nhosts){
      
      # loop across iterations
      for(j in 1:nrun){
        
        print(paste(i,j))
        sps.i <- sample(list.hosts,i) #random sample of i hosts
        subsdata <- subset(dat,HostCorrectedName%in%sps.i)
        
        # storing results
        rare.data[i,j] <- length(unique(subsdata$ParasiteCorrectedName))  
        rare.data[i,nrun+1] <- length(sps.i)/nhosts  
      }
    }
  }
  
  # method to count parasites within a subset of geographic coordinates
  if(by.host=='coords'){
    
    # generate array to store accumulation curves 
    rare.data <- array(NA,dim = c(ncoords,nrun+1))
    
    # loop along coordinates (i.e. locations)
    for(i in 1:ncoords){
      
      # loop along iterations
      for(j in 1:nrun){
        
        print(paste(i,j))
        coord.i <- sample(list.coords,i) #random sample of i coords
        position.coord <- which(list.coords%in%coord.i)
        subsdata <- dat[position.coord,]
        
        # storing results
        rare.data[i,j] <- length(unique(subsdata$ParasiteCorrectedName))  
        rare.data[i,nrun+1] <- length(coord.i)/ncoords
      }
    }
    
  }
  
  # method to count parasites within a subset of ecoregions
  if(by.host == 'ecoreg'){
    
    # generate array to store accumulation curves
    rare.data <- array(NA,dim = c(necoregs,nrun+1))
    
    # loop along ecoregions
    for(i in 1:necoregs){
      
      # loop along iterations
      for(j in 1:nrun){
        print(paste(i,j))
        ecoreg.i <- sample(list.ecoregs,i) #random sample of i ecoregions
        position.ecoreg <- which(list.ecoregs%in%ecoreg.i)
        subsdata <- dat[position.ecoreg,]
        
        # storing results
        rare.data[i,j] <- length(unique(subsdata$ParasiteCorrectedName))  
        rare.data[i,nrun+1] <- length(ecoreg.i)/necoregs
      }
    }
    
  }
  
  return(rare.data)
  
}


