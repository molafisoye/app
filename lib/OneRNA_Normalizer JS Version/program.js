
const csvParse = require('csv-parse');
const fs = require('fs');

(async (filename) => {
    const parse = csvParse({});
    const genes = {};
    const data = [];
    let values = [];

    parse.on('data', (breaker) => {
        data.push(breaker);
    });

    parse.on('end', async () => {
        const nCols = data.reduce((greater, { length }) => length > greater ? length : greater, 0);
       
        for (let col = 3; col < nCols; col++) {
            values.push([]);
            data.forEach((breaker, x) => {
                if (x === 0) {
                    return;
                }
                let str = breaker[col];

                const gene = breaker.slice(0, 3).join(',');
                
                const num = evenRound(parseFloat(str));
                str = formatFloat(num);

                if (!genes[gene]) {
                    genes[gene] = [num];
                    values[values.length-1].push(str);
                } else if(genes[gene].length < values.length) {
                    genes[gene].push(num);
                    values[values.length-1].push(str);
                }
            });
        }
        
        values.forEach(L => {
            L.sort();
        });
        
        console.log(data[0].join(','));
        
        let ite = 0;

        let entries = Object.entries(genes);
        entries.sort((a, b) => (
            a[0].localeCompare(b[0])
        ));
        
        for (let [value, L] of entries) {
            const TheCounts = [];
            for (let Y = 0; Y < L.length; Y++) {
                TheCounts[Y] = parseFloat(getNormalized(Y,parseInt(L[Y]), values, ite));
            }
            const CntStr = TheCounts.join(',');
            console.log(value+','+CntStr);
        }
    })

    fs.createReadStream(filename).pipe(parse);
})(process.argv[2]);

const getNormalized = (col, n, values, ite) => {
    const find = formatFloat(n);
    const i = values[col].indexOf(find);
    let sum = 0;

    if (i === -1) {
        // console.error(`Didn't find ${n} in col ${col}`);
    } else {
        for (let x=0; x < values.length; x++) {
            // if (ite === 1) {
            //     console.log('Avg', sum, values.length, sum/values.length)
            // }
            sum += parseInt(values[x][i]);
        }
    }


    return sum/values.length;
}

const formatFloat = (num) => String(Math.round(num)).padStart(10, '0');

const evenRound = (num, decimalPlaces = 0) => {
    const m = Math.pow(10, decimalPlaces);
    const n = +(decimalPlaces ? num * m : num).toFixed(12); // Avoid rounding errors
    const i = Math.floor(n), f = n - i;
    const e = 1e-8; // Allow for rounding errors in f
    const r = (f > 0.5 - e && f < 0.5 + e) ?
                ((i % 2 == 0) ? i : i + 1) : Math.round(n);
    return decimalPlaces ? r / m : r;
}