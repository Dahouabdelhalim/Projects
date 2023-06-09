# Script R para coletar os tweets sobre voto impresso entre 27/06/21 e 17/07/21

library(academictwitteR) # Carregando a biblioteca
# Mais informações sobre a biblioteca, acesse: https://CRAN.R-project.org/package=academictwitteR

# Para fazer o download de todos os tweets (n=3412)
tweets <- get_all_tweets(
    query = "voto impresso",
    start_tweets = "2021-06-27T00:00:00Z",
    end_tweets = "2021-07-17T23:59:59Z",
    n = Inf,
    data_path = "data2/",
    bind_tweets = FALSE,
    country= "br"
)

# Para fazer o download apenas de tweets, sem replies, retweets e quotes
onlytweets <-get_all_tweets(
    query = "voto impresso",
    start_tweets = "2021-06-27T00:00:00Z",
    end_tweets = "2021-07-17T23:59:59Z",
    n = Inf,
    data_path = "data2/",
    bind_tweets = FALSE,
    country= "br",
    is_retweet = FALSE,
    is_reply = FALSE,
    is_quote = FALSE
)

# Para converter os dados em uma base mais organizada visualmente
onlytweets <- bind_tweets(data_path = "data2/", output_format = "tidy")

# Fim