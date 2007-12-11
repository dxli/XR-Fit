// Copyright (c) 1999-2001 David Muse
// See the COPYING file for more information

#ifndef HTTP_H
#define HTTP_H

#include <inputoutput.h>

// The http class provides functions for generating http headers.

// http 1.1

class   http : public inputoutput {
        public:

                        http(void *apistruct);

                // the following enumeration types codes
                // for use in the functions below and cause the function
                // to send the corresponding code

                // enums
                typedef enum {
                        // success codes
                        ok,
                        created,
                        accepted,
                        partialinformation,
                        noresponse,
                        // error codes
                        badrequest,
                        unauthorized,
                        paymentrequired,
                        forbidden,
                        notfound,
                        internalerror,
                        notimplemented,
                        servicetemporarilyoverloaded,
                        gatewaytimeout,
                        // redirection codes
                        moved,
                        found,
                        method,
                        notmodified
                } statuscode;
                
                typedef enum {
                        get_request,
                        head_request,
                        post_request
                } requestmethod;
                
                typedef enum {
                        xgzip,
                        xcompress
                } encodingmethod;

                // The following functions are convenience functions for
                // sending http headers.  The name of the function is the
                // header that is sent and the arguments are for the
                // possible fields of the header.  A null or empty character
                // string will leave that field blank.
                
                // status line
                void    status(strstream *container,
                                char *serverprotocol, statuscode code);
                void    status(strstream *container,
                                char *protocol, char *protocolversion,
                                statuscode code);

                // general header fields
                void    date(strstream *container, char *httpdate);
                void    pragma(strstream *container, char *pragmadirective);
                void    pragma(strstream *container, char *token, char *word);

                // entity header fields
                void    allow(strstream *container, requestmethod method, ...);
                void    contentEncoding(strstream *container,
                                        encodingmethod method);
                void    contentLength(strstream *container, int length);
                void    contentType(strstream *container,
                                        char *mimetype, char *mimesubtype, 
                                        char *charset);
                void    multiPartContentType(strstream *container,
                                        char *subtype, char *charset, 
                                        char *boundary);
                void    expires(strstream *container,
                                        char *httpdate);
                void    lastModified(strstream *container,
                                        char *httpdate);
                

                // response header fields
                void    location(strstream *container, char *url);
                void    location(strstream *container, char *path, char *page);
                void    location(strstream *container,
                                char *protocol, char *host, char *port,
                                char *url);
                void    location(strstream *container,
                                char *protocol, char *host, char *port,
                                char *path, char *page);
                void    server(strstream *container,
                                char *product, char *version);
                void    wwwAuthenticate(strstream *container,
                                char *challenge, ...);
                void    setCookie(strstream *container,
                                char *name, char *value, char *path,
                                char *domain, char *expires, int secure);

                // request header fields
                void    authorization(strstream *container, char *credentials);
                void    from(strstream *container, char *mailbox);
                void    ifModifiedSince(strstream *container, char *httpdate);
                void    referer(strstream *container,
                                char *protocol, char *host, char *port,
                                char *path);
                void    userAgent(strstream *container,
                                char *product, char *version);


        private:
                #include <private/http.h>

};

#endif

