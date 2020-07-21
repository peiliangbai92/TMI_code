// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>



using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
mat soft_cpp(mat L, vec weight, double lambda){
    int width=L.n_cols;
    for(int x=0; x < width; x++) {
        double lambda_w = lambda*(1+ weight(x));
        if (L(0, x) > lambda_w)
            L(0, x) = L(0, x) - lambda_w;
        else if (L(0, x) <  - lambda_w)
            L(0, x) = L(0, x) + lambda_w;
        else
            L(0,x) = 0.0;
    }
    return L;
}

//[[Rcpp::export]]
mat pred_cpp(mat Y, mat phi, int p, int T, int k, int h){
    mat concat_Y (k,p+h, fill::zeros); 
    concat_Y(span::all , span(0, p-1)) = Y(span::all , span( T-p, T-1) );
    for( int j = 0; j < h; j++){
        mat temp (k, 1,fill::zeros );
        for(int i = 0; i < p; i ++){
            temp = temp + phi(span::all ,  span( i*k, (i+1)*k -1 )) * concat_Y(span::all , p+j-i-1);
        }
        concat_Y(span::all, p + j) = temp;
    }
    return concat_Y(span::all , p + h -1 );

}



// [[Rcpp::export]]
List var_lasso_brk(NumericMatrix data, NumericVector lambda, int p, int max_iteration, double tol ){
    int k = data.ncol(); int T = data.nrow(); int T_1 = T;
    int lambda_len = lambda.length();

    mat iter (k, lambda_len, fill::zeros); 
    mat phi_hat (k,k*p, fill::zeros); 
    mat phi_hat_fista (max_iteration, k*p, fill::zeros); 
    mat phi_hat_temp (k, k*p*lambda_len, fill::zeros); 
    vec pred_error =  zeros<vec>(lambda_len);

    mat data_m(data.begin(), T, k);
    mat Y = data_m.t();
    Y = Y(span::all, span(p, T_1-1));
    mat Z(k*p,T_1-p, fill::zeros);
    for(int i = 0; i < T_1 - p; i++) {
        for(int j = 1; j <= p; j++){
            Z(  span( (j-1)*k, j*k -1 ), i )= data_m( i+p-j, span::all).t();
        }

    }
    //Rcout << Z.n_rows << '\n'; 
    //Rcout << Z.n_cols << '\n'; 

    vec s = svd( Z );
    double step_size =1/pow(max(s) ,2); 

    vec weight = zeros<vec>(k*p);
    mat forecast_full (k, T_1-p, fill::zeros); 
    for (int ll = 0; ll < lambda_len; ll++ ){        
        for (int ii = 0; ii < k; ii ++ ){
            int l = 2;
            while( l < max_iteration){
                l  = l + 1;
                // Rcout << phi_hat_fista(l-2,span::all);
                mat phi_temp = phi_hat_fista(l-2,span::all) + ((l-2.0)/(l+1.0)) * (phi_hat_fista(l-2,span::all)  - phi_hat_fista(l-3,span::all) );
                //Rcout << phi_temp;
                mat phi_new =  phi_temp + step_size *  (Z  * (   Y( ii, span::all) - phi_temp*Z ).t() ).t();

                phi_new  =  soft_cpp(phi_new, weight , lambda(ll) );
                mat abs_temp = abs(phi_new - phi_temp);
                double max_temp = abs_temp.max();
                //Rprintf( "%f \n", max_temp);
                if ( max_temp < tol) {
                    phi_hat_temp( ii, span(  ll*k*p, (ll+1)*k*p -1 )) = phi_new;
                    break;
                } 
                if ( max_temp > tol ) {
                    // Rprintf( "%i \n", l);
                    // Rprintf( "%f \n", max_temp);
                    phi_hat_fista(l-1,span::all) =  phi_new;  
                    // Rcout << phi_new;
                }
            }           
            iter(ii,ll)  =  l;
            
        }

        //Rcout << phi_hat_temp;
        //Rcout << phi_hat_temp.n_cols;
        //Rcout << phi_hat_temp.n_rows;

        mat forecast (k, T_1-p, fill::zeros); 
        for(int j = p-1; j < T_1 - 1; j++){
        	//Rcout <<  pred_cpp( data_m.t(), phi_hat_temp(span::all ,  span( ll*k*p,(ll+1)*k*p-1) ) , p , j+1, k, 1); 
        	//Rcout << j; 
            forecast(span::all, (j-p+1)) = pred_cpp( data_m.t(), phi_hat_temp(span::all ,  span( ll*k*p,(ll+1)*k*p-1) ) , p , j+1, k, 1); 
            //Rcout <<  forecast(span::all, j); 
                     
        }
        
        // Rcout <<  data_m( span(p ,T_1-1), span::all);
        //Rcout <<  pow(  data_m( span(p ,T_1-1), span::all).t() - forecast , 2) ;
        forecast_full =forecast;
        pred_error(ll) = accu( pow(  (data_m( span(p ,T_1-1), span::all).t() - forecast) , 2) );
        //pred_error(ll) = accu( pow(  (data_m( span(p ,T_1-1), span::all).t() - forecast)/ (data_m( span(p ,T_1-1), span::all).t()), 2) );
        //Rcout << pred_error(ll);
        //Rcout << forecast;
        //Rcout<< pow(  (data_m( span(p ,T_1-1), span::all).t() - forecast)/ (data_m( span(p ,T_1-1), span::all).t()), 2);
    }
    int ll_final = 1;
    phi_hat  = phi_hat_temp(span::all , span( (ll_final-1)*k*p, ll_final*k*p -1 ));
    

    
    return List::create(Named("phi.hat")= phi_hat,Named("iter")= iter, Named("pred.error")= pred_error, Named("data_m")= data_m( span(p,T_1-1), span::all).t(), Named("forecast_full")=  forecast_full );
	
}


